import numpy as np
import xarray as xr
import os
import eddytools as et  # type: ignore
import pickle
import cftime
import datetime
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.simplefilter('ignore', UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

def detection_save_m(depth_index, yi, yf, Config, mp_cpu, Model, redo = True, old_interp = True):
    domain_path = Config['domain_path']
    prefixe_interpolation_c = Config['prefixe_interpolation_c']
    prefixe_interpolation = Config['prefixe_interpolation']
    prefixe_detection_m = Config['prefixe_detection_m']
    suffixe_all = Config['suffixe_all']
    suffixe_all_old = 'REF12_0'
    save_folder_path = Config['save_folder_path']
    save_folder_path_old_interp = '/home/nplanat/scratch/REF12/'

    ds_zgr = xr.open_dataset(domain_path)
    L_depth_value = ds_zgr.nav_lev.values.astype(int)
    depth = L_depth_value[depth_index]
    detection_parameters = Config['detection_parameters']
    if not os.path.exists(Config['datapath_c'] + prefixe_interpolation_c +'_depthi' + str(depth_index)+Config['detection_parameters']['suffixe_OW']+'.nc') :
        print('I cannot find mean OW_std, I look here :')
        print(Config['datapath_c'] + prefixe_interpolation_c +'_depthi' + str(depth_index)+Config['detection_parameters']['suffixe_OW']+'.nc')
    else:
        data_OW_std = xr.open_dataset(Config['datapath_c'] + prefixe_interpolation_c +'_depthi' + str(depth_index)+Config['detection_parameters']['suffixe_OW']+'.nc')
        if old_interp and not np.all(np.array([os.path.exists(save_folder_path_old_interp + prefixe_interpolation + str(year) +'_'+'{0:02}'.format(mo)+ '_depthi' + str(depth_index) + suffixe_all_old+'.nc') for year in np.arange(yi, yf+1) for mo in range(1, 13)])):
            print('I cannot find interpolation for that depth for all months/years, I look here : ')
            print(save_folder_path_old_interp + prefixe_interpolation + str(1995) +'_'+'{0:02}'.format(1)+ '_depthi' + str(depth_index) + suffixe_all_old+'.nc') 
        elif not old_interp and not np.all(np.array([os.path.exists(save_folder_path + prefixe_interpolation + str(year) +'_'+'{0:02}'.format(mo)+ '_depthi' + str(depth_index) + suffixe_all+'.nc') for year in np.arange(yi, yf+1) for mo in range(1, 13)])):
            print('I cannot find interpolation for that depth for all months/years, I look here : ')
            print(save_folder_path + prefixe_interpolation + str(1995) +'_'+'{0:02}'.format(1)+ '_depthi' + str(depth_index) + suffixe_all+'.nc') 
        else:
            for y in np.arange(yi, yf+1):
                year = str(y)   
                for mi in range(1, 13):
                    print(year, '/', mi)
                    first_day = str(cftime.datetime(int(year), mi, 1, calendar = detection_parameters['calendar']))[0:10]
                    
                    detection_parameters['start_time'] =first_day, # time range start
                    detection_parameters['lon1'] = Model['lon1']
                    detection_parameters['lon2'] = Model['lon2']
                    if mi==12:
                        last_day = year+'-12-31'
                    else:
                        last_day = str(cftime.datetime(y, mi+1, 1, calendar = detection_parameters['calendar'])-datetime.timedelta(days = 1))[0:10]
                    detection_parameters['end_time'] = last_day, # time range end
                    detection_parameters['start_time'] = detection_parameters['start_time'][0]
                    detection_parameters['end_time'] = detection_parameters['end_time'][0]
                    detection_parameters['min_dep']= depth-1, # minimum ocean depth where to look for eddies in m
                    detection_parameters['min_dep'] = detection_parameters['min_dep'][0]
                    if not os.path.exists(save_folder_path + prefixe_detection_m +  year+'_'+'{0:02}'.format(mi) +'_depthi' + str(depth_index) + suffixe_all+'_eddies.pickle'):
                        print(year)
                        #branching from REF12_2 : 
                        if old_interp:
                            save_datapath = save_folder_path_old_interp + prefixe_interpolation + year +'_'+'{0:02}'.format(mi)+ '_depthi' + str(depth_index) + suffixe_all_old+'.nc'
                        #normaly use this : 
                        else:
                            save_datapath = save_folder_path + prefixe_interpolation + year+'_'+'{0:02}'.format(mi)+ '_depthi' + str(depth_index) + suffixe_all+'.nc'
                        data_int = xr.open_mfdataset(save_datapath)
                        data_int = data_int.drop({'OW_std'})
                        data_int['OW_std'] = data_OW_std.OW_std.isel(time = 0)
                        detection_parameters['OW_thr'] = data_OW_std['OW_std'].values,
                        eddies = et.detection.detect_OW(data_int, detection_parameters, 'OW', 'vort', regrid_avoided=True, use_mp=True, mp_cpu=mp_cpu) 
                        with open(save_folder_path + prefixe_detection_m +  year+'_'+'{0:02}'.format(mi)+ '_depthi' + str(depth_index) + suffixe_all+'_eddies.pickle', 'wb') as f:
                            pickle.dump(eddies, f, pickle.HIGHEST_PROTOCOL)
                        f.close()
                    else:
                        print(year, '/', mi, 'depth: ', depth_index,  'I already have the detection for that depth ')
    return None
