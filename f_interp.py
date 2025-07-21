from xorca.lib import load_xorca_dataset # type: ignore
import dask # type: ignore
import xarray as xr
import eddytools as et # type: ignore
import numpy as np
from glob import glob
import f_recreate_CREG
import os
import datetime
import cftime
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.simplefilter('ignore', UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
from scipy.signal import convolve  

def nanconv(a, k, MODE):
    # Flat function for comparison.
    o = np.ones(np.shape(a))
    # Flat function with NaNs for comparison.
    on = np.ones(np.shape(a))
    # Find all the NaNs in the input and replace with 0
    an = np.where(~np.isnan(a), a, 0)
    on = np.where(~np.isnan(a), on, 0)
    # Calculate what a 'flat' function looks like after convolution.
    flat = convolve(on, k, mode=MODE)
    #
    # The line above will automatically include a correction for edge 
    # effects,
    # so remove that correction if the user does not want it.
    flat = flat / convolve(o, k, mode=MODE)
    #
    # Do the actual convolution
    output = convolve(an, k, mode=MODE) / flat
    return np.where(~np.isnan(a), output, np.nan)  

def spatial_std(data, wx, wy):
    window = np.ones((1, wy, wx)) / (wy * wx)
    ext = np.zeros((np.shape(data)[0], 
                    np.shape(data)[1] + wy, np.shape(data)[2] + wx))
    ext[:, int(wy/2):-int(wy/2), int(wx/2):-int(wx/2)] = data
    ext[:, 0:int(wy/2), :] = ext[:, int(wy):int(wy/2):-1, :]
    ext[:, -int(wy/2)::, :] = ext[:, -int(wy/2):-int(wy):-1, :]
    ext[:, :, 0:int(wx/2)] = ext[:, :, int(wx):int(wx/2):-1]
    ext[:, :, -int(wx/2)::] = ext[:, :, -int(wx/2):-int(wx):-1]
    std_tmp1 = np.abs(ext - nanconv(ext, window, "same")) ** 2
    std_tmp2 = nanconv(std_tmp1, window, "same") ** 0.5
    output =  xr.DataArray(std_tmp2[:, int(wy/2):-int(wy/2), int(wx/2):-int(wx/2)], 
                           coords=data.coords, dims=data.dims).mean("time") 
    return output  

def interpolate_new_rolling(depth_index, yi, yf, Config, Model, lsig = 50, cs_vert = 75, cs_x = 1580, cs_y = 1801, should_convert_cftime = False):
    domain_path = Config['domain_path']
    prefixe_interpolation = Config['prefixe_interpolation']
    suffixe_all = Config['suffixe_all']
    save_folder_path = Config['save_folder_path']
    datapath = Config['datapath']
    meshpath = Config['meshpath']
    maskpath = Config['maskpath']
    interpolation_parameters = Config['interpolation_parameters']
    ds_zgr = xr.open_dataset(domain_path)
    L_depth_value = ds_zgr.nav_lev.values.astype(int)
    depth = L_depth_value[depth_index]

    for y in np.arange(yi, yf+1):
        year = str(y)
        for mo in range(1, 13):
            #get lastday of month
            first_day = str(cftime.datetime(int(year), mo, 1, calendar = interpolation_parameters['calendar']))[0:10]
            if mo==12:
                last_day = year+'-12-31'
            else:
                last_day = str(cftime.datetime(int(year), mo+1, 1, calendar = interpolation_parameters['calendar'])-datetime.timedelta(days = 1))[0:10]

            interpolation_parameters['start_time'] = first_day, # time range start
            interpolation_parameters['end_time'] = last_day, # time range end
            interpolation_parameters['end_time'] = interpolation_parameters['end_time'][0]
            interpolation_parameters['start_time'] = interpolation_parameters['start_time'][0]
            if not os.path.exists(save_folder_path + prefixe_interpolation + year +'_'+'{0:02}'.format(mo)+ '_depthi' + str(depth_index) + suffixe_all+'.nc'):
                print('I do ', year, '/', str(mo))
                data_in = sorted(glob(datapath +Config['prefixe_data']+year+'m'+'*.nc'))
                print(data_in)
                with dask.config.set(**{'array.slicing.split_large_chunks': False}):
                    data = load_xorca_dataset(data_files=data_in, aux_files=[meshpath, maskpath], model_config='NEST',decode_cf = True,
                                            input_ds_chunks = {"time_counter": 1, "t": 1,
                                                                "z": cs_vert, "deptht": cs_vert, "depthu": cs_vert, "depthv": cs_vert, "depthw": cs_vert,
                                                                "x": cs_x, "y": cs_y},
                                            target_ds_chunks = {"t":30, "z_c": 1, "z_l": 1,
                                                                "x_c": cs_x, "x_r": cs_x, "y_c": cs_y, "y_r": cs_y})
                print(data.dims)
                if should_convert_cftime:
                    Times = []
                    for i in range(len(data.t)):
                        Times.append(cftime.datetime(data.t.dt.year[i], data.t.dt.month[i], data.t.dt.day[i], data.t.dt.hour[i], calendar = Model['calendar']))
                    data['t'] = ('t', Times)

                data_z = data.isel(z_c=depth_index, z_l=depth_index)
                data_y, grid, metrics = f_recreate_CREG.do_all(data_z.sel(t=slice( interpolation_parameters['start_time'], interpolation_parameters['end_time'])), domain_path)

                # Calculate vorticity and Okubo-Weiss parameter and make sure the chunk sizes are as before.
                data_OW = et.okuboweiss.calc(data_y, grid, 'vozocrtx', 'vomecrty').chunk({'x_c': cs_x, 'x_r': cs_x,'y_c': cs_y, 'y_r': cs_y})
                # Merge the new variables `OW` and `vort` to the dataset `data`
                data_y = xr.merge([data_y, data_OW], compat='override')
                # get rid of halo-points
                data_y = data_y.isel(y_c=slice(0, -1), y_r=slice(0, -1), x_c=slice(None, -2), x_r=slice(None, -2))
                # extract region really needed for analysis
                data_cut = data_y.isel(x_c=slice(Model['i1'], Model['i2']), x_r=slice(Model['i1'], Model['i2']),
                                    y_c=slice(Model['j1'], Model['j2']), y_r=slice(Model['j1'], Model['j2']))
                # Define the parameters for the interpolation
                avoid_regrid = True
                data_int, regridder = et.interp.horizontal(data_cut, metrics, interpolation_parameters, weights=None, avoid_regrid=avoid_regrid)

                OW_tmp = data_int['OW']#.sel(time=slice(time_range_start, time_range_end))
                OW_tmp = OW_tmp.where(OW_tmp != 0).chunk({"time": 1}).persist()

                if avoid_regrid == True:
                    x = 'x'
                    y = 'y'
                else:
                    x = 'lon'
                    y = 'lat'

                lon_tmp = OW_tmp[x].where(OW_tmp[x] > 0, other=OW_tmp[x] + 360.)
                OW_tmp = OW_tmp.assign_coords({x: lon_tmp})

                #data_int = data_int.drop({'OW_std'})
                data_int['OW_std'] = spatial_std(OW_tmp, lsig, lsig).where(data_int.tmask > 0)
                #data_int['OW_std'] = ([y, x], mean_OW_spatial_std.where(data_int.tmask > 0)))
                #data_int["OW_std"] = data_int["OW_std"].where(data_int.tmask > 0)
                save_datapath= save_folder_path + prefixe_interpolation  + year +'_'+'{0:02}'.format(mo)+ '_depthi' + str(depth_index) + suffixe_all+'.nc'
                data_int.to_netcdf(save_datapath)                    
            else:
                print(year, '/', str(mo), ': Interp done ')    
    return None
