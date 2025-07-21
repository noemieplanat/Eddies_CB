from glob import glob
import xarray as xr
import os
import numpy as np

def save_grid_REF12(Config):
    prefixe_interpolation = Config['prefixe_interpolation']
    suffixe_all = Config['suffixe_all']
    save_folder_path = Config['save_folder_path']

    if not os.path.exists(save_folder_path+'Grid_selection_Area.nc'):
        path_interp= sorted(glob(save_folder_path + prefixe_interpolation +'*'+ suffixe_all+'.nc'))
        print(path_interp)
        data_int = xr.open_mfdataset(path_interp[0])
        ds_file = data_int.drop({'OW', 'vort', 'OW_std'}).isel(time = 0)
        ds_file.to_netcdf(save_folder_path+'Grid_selection_Area.nc')
    ds = xr.open_dataset(save_folder_path+'Grid_selection_Area.nc')
    return ds

def do_empty_array_REF12_fpoints(Config):
    ''' conserves Halo points if they exist'''
    save_folder_path = Config['save_folder_path']
    domain_path = Config['domain_path']
    if not os.path.exists(save_folder_path + 'Empty_f.nc'):
        de = xr.open_dataset(domain_path)
        de['x'] = de.x
        de['y'] = de.y
        de['z'] = de.z
        de = de.drop({'bathy_meter', 'bottom_level', 'e1f', 'e2f', 'e1u', 'e2u', 'e1t', 'e2t', 'e1v', 'e2v', 'e3f_0', 'e3t_0', \
                    'e3t_1d', 'e3u_0', 'e3v_0', 'e3uw_0','e3vw_0', 'e3w_0', 'e3w_1d', 'ff_f', 'ff_t', 'glamt', 'glamu',\
                    'glamv', 'gphiu', 'gphiv', 'gphit', 'nav_lat', 'nav_lon', 'jperio', 'jpiglo', 'jpjglo', 'jpkglo', 'ln_isfcav', 'ln_sco',\
                    'ln_zco','ln_zps', 'namelist_cfg', 'top_level'})
        de = de.rename({'glamf':'nav_lon', 'gphif':'nav_lat'})
        de.to_netcdf(save_folder_path+'Empty_f.nc')  
    de = xr.open_dataset(save_folder_path+ 'Empty_f.nc')
    return de


def get_distance_latlon2m(lon1, lat1, lon2, lat2):
    '''Get the distance in km between two points (lon1, lat1) and (lon2, lat2)
    based on longitude and latitude.
    We use the haversine approach which assumes the Earth is a sphere and can
    induce errors of up to 0.5%. (That means if we look for eddies in a search
    radius of 50 km, the actual radius could be between 49.75 km and 50.25 km.
    Errors due to interpolation etc are probably much larger!)

    Parameters
    ----------
    lon1 : float or int
        Longitude of the first point.
    lat1 : float or int
        Latitude of the first point.
    lon1 : float or int
        Longitude of the first point.
    lat1 : float or int
        Latitude of the first point.

    Returns
    -------
    d : float or int
        Distance between the two points (lon1, lat1) and (lon2, lat2)
    '''
    # Approximate radius of the Earth in km
    R = 6371.
    lon1 = np.radians(lon1)
    lat1 = np.radians(lat1)
    lon2 = np.radians(lon2)
    lat2 = np.radians(lat2)
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = ((np.sin(dlat / 2.) ** 2.) + np.cos(lat1)
         * np.cos(lat2) * (np.sin(dlon / 2.) ** 2.))
    c = 2 * np.arcsin(np.sqrt(a))
    return R * c * 1e3