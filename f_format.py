import os
import numpy as np
import xarray as xr
import pickle

def format_save_tracked_filtered(Config, ind_depth, filter):
    print(Config['save_folder_path'])
    print(Config['tracking_formated_prefixe'])
    print(Config['filtering_parameters'][filter]['name'])
    if not os.path.exists(Config['save_folder_path']+Config['tracking_formated_prefixe']+'_'+Config['filtering_parameters'][filter]['name']+'_'+str(ind_depth)+'_'+Config['suffixe_all']+'.nc'):
        if os.path.exists(Config['save_folder_path'] +Config['prefixe_tracking']+str(ind_depth)+Config['suffixe_all']+'_tracking.pickle') and os.path.exists(Config['save_folder_path'] +Config['prefixe_tracking']+str(ind_depth)+Config['suffixe_all']+'_tracking_1d.pickle'):

            tracked = load_tracked(Config['save_folder_path'], ind_depth, Config['prefixe_tracking'], Config['suffixe_all'])
            tracked_1d = load_tracked_1d(Config['save_folder_path'], ind_depth, Config['prefixe_tracking'], Config['suffixe_all'])
            
            formated = format_xarray_tracked(tracked, tracked_1d)
            formated = formated.where(np.logical_or(formated.lon_i<=Config['filtering_parameters'][filter]['lon1'], formated.lon_i>=Config['filtering_parameters'][filter]['lon2']), drop = True)
            formated = formated.where(formated.lat_i>=Config['filtering_parameters'][filter]['lat1'], drop = True).where(formated.lat_i<=Config['filtering_parameters'][filter]['lat2'], drop = True)
            formated = formated.where(formated.depth_i>=Config['filtering_parameters'][filter]['bathy_min'], drop = True)
            formated = formated.where(np.logical_or(formated.lon_f<=Config['filtering_parameters'][filter]['lon1'], formated.lon_f>=Config['filtering_parameters'][filter]['lon2']), drop = True)
            formated = formated.where(formated.lat_f>=Config['filtering_parameters'][filter]['lat1'], drop = True).where(formated.lat_f<=Config['filtering_parameters'][filter]['lat2'], drop = True)
            formated = formated.where(formated.depth_f>=Config['filtering_parameters'][filter]['bathy_min'], drop = True)
            formated = formated.where(formated.duration>=Config['filtering_parameters'][filter]['duration'], drop = True)
            formated.to_netcdf(Config['save_folder_path']+Config['tracking_formated_prefixe']+'_'+Config['filtering_parameters'][filter]['name']+'_'+str(ind_depth)+'TEMP.nc')
            format = xr.open_dataset(Config['save_folder_path']+Config['tracking_formated_prefixe']+'_'+Config['filtering_parameters'][filter]['name']+'_'+str(ind_depth)+'TEMP.nc')
            if len(format)>1 and len(format.eddy_num.where(format.time_i.dt.year>=Config['filtering_parameters'][filter]['yi'], drop = True).where(format.time_f.dt.year<=Config['filtering_parameters'][filter]['yf'], drop = True))>1:
                format = format.where(format.time_i.dt.year>=Config['filtering_parameters'][filter]['yi'], drop = True).where(format.time_f.dt.year<=Config['filtering_parameters'][filter]['yf'], drop = True)
                format['type'] = format['type'].astype(str)
                print('I have formatted : ', format.dims['eddy_num'])
            else:
                print('len tracked not good, depthi', ind_depth)
            format.to_netcdf(Config['save_folder_path']+Config['tracking_formated_prefixe']+'_'+Config['filtering_parameters'][filter]['name']+'_'+str(ind_depth)+'_'+Config['suffixe_all']+'.nc')
            os.remove(Config['save_folder_path']+Config['tracking_formated_prefixe']+'_'+Config['filtering_parameters'][filter]['name']+'_'+str(ind_depth)+'TEMP.nc')

            print('ind_depth ', ind_depth, ' done')
        else:
            print('ind_depth ', ind_depth, 'missing track')
    else:
        print('ind_depth ', ind_depth, ' already done')
    return None
def format_xarray_tracked(tracked, tracked_1d):
    c = len(tracked)
    d = len(tracked_1d)
    print('I have tracked ', c+d, ' eddies at that depth')
    ds = xr.Dataset(coords = {'eddy_num':np.arange(c+d)})
    LON_i = []
    LAT_i = []
    TIME_i = []
    LON_f = []
    LAT_f = []
    TIME_f = []
    LON_m = []
    LAT_m = []
    DURATION = []
    AMP = []
    AMP_i = []
    AMP_f = []
    AREA = []
    AREA_i = []
    AREA_f = []
    TYPE = []
    SCALE = []
    DEPTH_i = []
    DEPTH_f = []
    for i in range(len(tracked)):
        #if i%int(len(tracked)/10)==0:
        #    print(np.round(i/len(tracked)*100, 1), '%')
        LON_i.append(float(tracked[i]['lon'][0]))
        LAT_i.append(float(tracked[i]['lat'][0]))
        TIME_i.append(tracked[i]['time'][0])
        LON_f.append(float(tracked[i]['lon'][-1]))
        LAT_f.append(float(tracked[i]['lat'][-1]))
        TIME_f.append(tracked[i]['time'][-1])
        DURATION.append(len(tracked[i]['time']))
        AMP.append(np.mean(tracked[i]['amp']))
        AMP_i.append(float(tracked[i]['amp'][0]))
        AMP_f.append(float(tracked[i]['amp'][-1]))
        AREA.append(np.mean(tracked[i]['area']))
        DEPTH_i.append(float(tracked[i]['depth'][0]))
        DEPTH_f.append(float(tracked[i]['depth'][-1]))
        AREA_i.append(float(tracked[i]['area'][0]))
        AREA_f.append(float(tracked[i]['area'][-1]))
        TYPE.append(tracked[i]['type'])
        SCALE.append(np.mean(tracked[i]['scale']))
        
    for i in range(len(tracked_1d)):
        #if i%int(len(tracked_1d)/10)==0:
        #    print(np.round(i/len(tracked_1d)*100, 1), '%')
        LON_i.append(float(tracked_1d[i]['lon']))
        LAT_i.append(float(tracked_1d[i]['lat']))
        TIME_i.append(tracked_1d[i]['time'].tolist())
        LON_f.append(float(tracked_1d[i]['lon']))
        LAT_f.append(float(tracked_1d[i]['lat']))
        TIME_f.append(tracked_1d[i]['time'])
        DURATION.append(1)
        AMP.append(np.mean(tracked_1d[i]['amp']))
        AMP_i.append(float(tracked_1d[i]['amp']))
        AMP_f.append(float(tracked_1d[i]['amp']))
        AREA.append(np.mean(tracked_1d[i]['area']))
        DEPTH_i.append(float(tracked_1d[i]['depth'].mean()))
        DEPTH_f.append(float(tracked_1d[i]['depth'].mean()))
        AREA_i.append(float(tracked_1d[i]['area']))
        AREA_f.append(float(tracked_1d[i]['area']))
        TYPE.append(tracked_1d[i]['type'])
        SCALE.append(np.mean(tracked_1d[i]['scale']))

    ds['lon_i'] = (('eddy_num'), np.array(LON_i))  
    ds['lat_i'] = (('eddy_num'), np.array(LAT_i))   
    ds['time_i'] = (('eddy_num'), TIME_i)
    ds['amp_m'] = (('eddy_num'), np.array(AMP))   
    ds['amp_i'] = (('eddy_num'), np.array(AMP_i)) 
    ds['amp_f'] = (('eddy_num'), np.array(AMP_f)) 
    ds['area_m'] = (('eddy_num'), np.array(AREA))  
    ds['area_i'] = (('eddy_num'), np.array(AREA_i)) 
    ds['area_f'] = (('eddy_num'), np.array(AREA_f)) 
    ds['type'] = (('eddy_num'), TYPE)
    ds['depth_i'] = (('eddy_num'), np.array(DEPTH_i))
    ds['depth_f'] = (('eddy_num'), np.array(DEPTH_f))
    ds['scale_m'] = (('eddy_num'), np.array(SCALE))   
    ds['duration'] = (('eddy_num'), np.array(DURATION))   
    ds['lon_f'] = (('eddy_num'), np.array(LON_f))   
    ds['lat_f'] = (('eddy_num'), np.array(LAT_f))   
    ds['time_f'] = (('eddy_num'), np.array(TIME_f))  
    return ds
def load_tracked_1d(path, depth_ind, prefixe_tracking, suffixe_all):
    '''returns a dictionary with one item per detected eddy'''
    file=path +prefixe_tracking+str(depth_ind)+suffixe_all+'_tracking_1d.pickle'
    infile = open(file,'rb')
    tracked = pickle.load(infile)
    infile.close()
    return tracked
def load_tracked(path, depth_ind, prefixe_tracking, suffixe_all):
    '''returns a dictionary with one item per detected eddy'''
    files=path +prefixe_tracking+str(depth_ind)+suffixe_all+'_tracking.pickle'
    infile = open(files,'rb')
    tracked = pickle.load(infile)
    infile.close()
    return tracked
def load_tracked_properties(path, depth_ind, prefixe_properties, suffixe_all):
    file = path + prefixe_properties+'tracked_depthi' + str(depth_ind) + suffixe_all+'tracked_properties.pickle'
    infile = open(file,'rb')
    tracked = pickle.load(infile)
    infile.close()
    print('prop', len(tracked))
    print('file', file)
    return tracked
def load_tracked_properties_1d(path, depth_ind, prefixe_properties, suffixe_all):
    file = path + prefixe_properties+'tracked_depthi' + str(depth_ind) + suffixe_all+'tracked_properties_1d.pickle'
    infile = open(file,'rb')
    tracked = pickle.load(infile)
    infile.close()
    print('prop_1d', len(tracked))
    print('file 1d', file)
    return tracked
def format_xarray_tracked_with_properties(tracked, tracked_1d, properties_tracked, properties_tracked_1d, Liste_properties_to_track, duration = 1):
    if duration >=2:
        c = len(tracked)
        if c != len(properties_tracked):
            print('I have a problem with tracked properties !!!!!---------------')
            print('tracked', len(tracked))
            print('tracked p ', len(properties_tracked))
            return None
        else:
            print('I have tracked ', c, ' eddies at that depth')
            ds = xr.Dataset(coords = {'eddy_num':np.arange(c)})
            LON_i = []
            LAT_i = []
            TIME_i = []
            LON_f = []
            LAT_f = []
            TIME_f = []
            LON_m = []
            LAT_m = []
            DURATION = []
            AMP = []
            AMP_i = []
            AMP_f = []
            AREA = []
            AREA_i = []
            AREA_f = []
            TYPE = []
            SCALE = []
            DEPTH_i = []
            DEPTH_f = []
            PROP_i = [[] for x in range(len(Liste_properties_to_track))]
            PROP_f = [[] for x in range(len(Liste_properties_to_track))]
            for i in range(len(properties_tracked)):
                #if i%int(len(tracked)/10)==0:
                #    print(np.round(i/len(tracked)*100, 1), '%')
                LON_i.append(float(tracked[i]['lon'][0]))
                LAT_i.append(float(tracked[i]['lat'][0]))
                TIME_i.append(tracked[i]['time'][0])
                LON_f.append(float(tracked[i]['lon'][-1]))
                LAT_f.append(float(tracked[i]['lat'][-1]))
                TIME_f.append(tracked[i]['time'][-1])
                DURATION.append(len(tracked[i]['time']))
                AMP.append(np.mean(tracked[i]['amp']))
                AMP_i.append(float(tracked[i]['amp'][0]))
                AMP_f.append(float(tracked[i]['amp'][-1]))
                AREA.append(np.mean(tracked[i]['area']))
                AREA_i.append(float(tracked[i]['area'][0]))
                DEPTH_i.append(float(tracked[i]['depth'][0]))
                DEPTH_f.append(float(tracked[i]['depth'][-1]))
                AREA_f.append(float(tracked[i]['area'][-1]))
                TYPE.append(tracked[i]['type'])
                SCALE.append(np.mean(tracked[i]['scale']))
                for (ip, p) in zip(np.arange(len(Liste_properties_to_track)), Liste_properties_to_track):
                    PROP_i[ip].append(float(properties_tracked[i][p][0]))
                    PROP_f[ip].append(float(properties_tracked[i][p][-1]))

            ds['lon_i'] = (('eddy_num'), np.array(LON_i))  
            ds['lat_i'] = (('eddy_num'), np.array(LAT_i))   
            ds['time_i'] = (('eddy_num'), TIME_i)
            ds['amp_m'] = (('eddy_num'), np.array(AMP))   
            ds['amp_i'] = (('eddy_num'), np.array(AMP_i)) 
            ds['amp_f'] = (('eddy_num'), np.array(AMP_f)) 
            ds['area_m'] = (('eddy_num'), np.array(AREA))  
            ds['area_i'] = (('eddy_num'), np.array(AREA_i)) 
            ds['area_f'] = (('eddy_num'), np.array(AREA_f)) 
            ds['type'] = (('eddy_num'), TYPE)
            ds['scale_m'] = (('eddy_num'), np.array(SCALE))   
            ds['duration'] = (('eddy_num'), np.array(DURATION))   
            ds['lon_f'] = (('eddy_num'), np.array(LON_f))   
            ds['lat_f'] = (('eddy_num'), np.array(LAT_f))   
            ds['time_f'] = (('eddy_num'), np.array(TIME_f))      
            ds['depth_i'] = (('eddy_num'), np.array(DEPTH_i))
            ds['depth_f'] = (('eddy_num'), np.array(DEPTH_f))
            for (ip, p) in zip(np.arange(len(Liste_properties_to_track)), Liste_properties_to_track):
                ds[Liste_properties_to_track[ip]+'_i'] = (('eddy_num'), np.array(PROP_i[ip]))
                ds[Liste_properties_to_track[ip]+'_f'] = (('eddy_num'), np.array(PROP_f[ip]))
            return ds
    else:
        c = len(tracked)+len(tracked_1d)
        if c != len(properties_tracked)+len(properties_tracked_1d):
            print('I have a problem with tracked properties !!!!!---------------')
            print('tracked', len(tracked))
            print('tracked 1d ', len(tracked_1d))
            print('tracked p ', len(properties_tracked))
            print('tracked p 1d ', len(properties_tracked_1d))
            return None
        else:
            print('I have tracked ', c, ' eddies at that depth')
            ds = xr.Dataset(coords = {'eddy_num':np.arange(c)})
            LON_i = []
            LAT_i = []
            TIME_i = []
            LON_f = []
            LAT_f = []
            TIME_f = []
            LON_m = []
            LAT_m = []
            DURATION = []
            AMP = []
            AMP_i = []
            AMP_f = []
            AREA = []
            AREA_i = []
            AREA_f = []
            TYPE = []
            SCALE = []
            DEPTH_i = []
            DEPTH_f = []
            PROP_i = [[] for x in range(len(Liste_properties_to_track))]
            PROP_f = [[] for x in range(len(Liste_properties_to_track))]
            for i in range(len(properties_tracked)):
                #if i%int(len(tracked)/10)==0:
                #    print(np.round(i/len(tracked)*100, 1), '%')
                LON_i.append(float(tracked[i]['lon'][0]))
                LAT_i.append(float(tracked[i]['lat'][0]))
                TIME_i.append(tracked[i]['time'][0])
                LON_f.append(float(tracked[i]['lon'][-1]))
                LAT_f.append(float(tracked[i]['lat'][-1]))
                TIME_f.append(tracked[i]['time'][-1])
                DURATION.append(len(tracked[i]['time']))
                AMP.append(np.mean(tracked[i]['amp']))
                AMP_i.append(float(tracked[i]['amp'][0]))
                AMP_f.append(float(tracked[i]['amp'][-1]))
                AREA.append(np.mean(tracked[i]['area']))
                AREA_i.append(float(tracked[i]['area'][0]))
                DEPTH_i.append(float(tracked[i]['depth'][0]))
                DEPTH_f.append(float(tracked[i]['depth'][-1]))
                AREA_f.append(float(tracked[i]['area'][-1]))
                TYPE.append(tracked[i]['type'])
                SCALE.append(np.mean(tracked[i]['scale']))
                for (ip, p) in zip(np.arange(len(Liste_properties_to_track)), Liste_properties_to_track):
                    PROP_i[ip].append(float(properties_tracked[i][p][0]))
                    PROP_f[ip].append(float(properties_tracked[i][p][-1]))


            for i in range(len(properties_tracked_1d)):
                #if i%int(len(tracked_1d)/10)==0:
                #    print(np.round(i/len(tracked_1d)*100, 1), '%')
                LON_i.append(float(tracked_1d[i]['lon']))
                LAT_i.append(float(tracked_1d[i]['lat']))
                TIME_i.append(tracked_1d[i]['time'].tolist())
                LON_f.append(float(tracked_1d[i]['lon']))
                LAT_f.append(float(tracked_1d[i]['lat']))
                TIME_f.append(tracked_1d[i]['time'])
                DURATION.append(1)
                AMP.append(np.mean(tracked_1d[i]['amp']))
                AMP_i.append(float(tracked_1d[i]['amp']))
                AMP_f.append(float(tracked_1d[i]['amp']))
                AREA.append(np.mean(tracked_1d[i]['area']))
                DEPTH_i.append(float(tracked_1d[i]['depth'].mean()))
                DEPTH_f.append(float(tracked_1d[i]['depth'].mean()))
                AREA_i.append(float(tracked_1d[i]['area']))
                AREA_f.append(float(tracked_1d[i]['area']))
                TYPE.append(tracked_1d[i]['type'])
                SCALE.append(np.mean(tracked_1d[i]['scale']))

                for (ip, p) in zip(np.arange(len(Liste_properties_to_track)), Liste_properties_to_track):
                    PROP_i[ip].append(float(properties_tracked_1d[i][p]))
                    PROP_f[ip].append(float(properties_tracked_1d[i][p]))
            ds['lon_i'] = (('eddy_num'), np.array(LON_i))  
            ds['lat_i'] = (('eddy_num'), np.array(LAT_i))   
            ds['time_i'] = (('eddy_num'), TIME_i)
            ds['amp_m'] = (('eddy_num'), np.array(AMP))   
            ds['amp_i'] = (('eddy_num'), np.array(AMP_i)) 
            ds['amp_f'] = (('eddy_num'), np.array(AMP_f)) 
            ds['area_m'] = (('eddy_num'), np.array(AREA))  
            ds['area_i'] = (('eddy_num'), np.array(AREA_i)) 
            ds['area_f'] = (('eddy_num'), np.array(AREA_f)) 
            ds['type'] = (('eddy_num'), TYPE)
            ds['scale_m'] = (('eddy_num'), np.array(SCALE))   
            ds['duration'] = (('eddy_num'), np.array(DURATION))   
            ds['lon_f'] = (('eddy_num'), np.array(LON_f))   
            ds['lat_f'] = (('eddy_num'), np.array(LAT_f))   
            ds['time_f'] = (('eddy_num'), np.array(TIME_f))      
            ds['depth_i'] = (('eddy_num'), np.array(DEPTH_i))
            ds['depth_f'] = (('eddy_num'), np.array(DEPTH_f))
            for (ip, p) in zip(np.arange(len(Liste_properties_to_track)), Liste_properties_to_track):
                ds[Liste_properties_to_track[ip]+'_i'] = (('eddy_num'), np.array(PROP_i[ip]))
                ds[Liste_properties_to_track[ip]+'_f'] = (('eddy_num'), np.array(PROP_f[ip]))
            return ds
def format_save_tracked_prop_filtered(Config, ind_depth, filter):
    if not os.path.exists(Config['save_folder_path']+Config['tracking_formated_prefixe']+'Prop_'+Config['filtering_parameters'][filter]['name']+'_'+str(ind_depth)+'_'+Config['suffixe_all']+'.nc'):
        if os.path.exists(Config['save_folder_path'] + Config['prefixe_properties_tracking']+'tracked_depthi' + str(ind_depth) + Config['suffixe_all']+'tracked_properties.pickle'):
            if Config['filtering_parameters'][filter]['duration']>1:
                tracked = load_tracked(Config['save_folder_path'], ind_depth, Config['prefixe_tracking'], Config['suffixe_all'])
                prop = load_tracked_properties(Config['save_folder_path'], ind_depth, Config['prefixe_properties_tracking'], Config['suffixe_all'])
                formated = format_xarray_tracked_with_properties(tracked, None, prop, None, Config['properties_parameters']['Properties_extract'], duration = Config['filtering_parameters'][filter]['duration'])
            else:
                prop_1d = load_tracked_properties_1d(Config['save_folder_path'], ind_depth, Config['prefixe_properties_tracking'], Config['suffixe_all'])
                tracked_1d = load_tracked_1d(Config['save_folder_path'], ind_depth, Config['prefixe_tracking'], Config['suffixe_all'])
                tracked = load_tracked(Config['save_folder_path'], ind_depth, Config['prefixe_tracking'], Config['suffixe_all'])
                prop = load_tracked_properties(Config['save_folder_path'], ind_depth, Config['prefixe_properties_tracking'], Config['suffixe_all'])
                formated = format_xarray_tracked_with_properties(tracked, tracked_1d, prop, prop_1d, Config['properties_parameters']['Properties_extract'], duration = Config['filtering_parameters'][filter]['duration'])
            
            formated = formated.where(np.logical_or(formated.lon_i<=Config['filtering_parameters'][filter]['lon1'], formated.lon_i>=Config['filtering_parameters'][filter]['lon2']), drop = True)
            formated = formated.where(formated.lat_i>=Config['filtering_parameters'][filter]['lat1'], drop = True).where(formated.lat_i<=Config['filtering_parameters'][filter]['lat2'], drop = True)
            formated = formated.where(formated.depth_i>=Config['filtering_parameters'][filter]['bathy_min'], drop = True)
            formated = formated.where(np.logical_or(formated.lon_f<=Config['filtering_parameters'][filter]['lon1'], formated.lon_f>=Config['filtering_parameters'][filter]['lon2']), drop = True)
            formated = formated.where(formated.lat_f>=Config['filtering_parameters'][filter]['lat1'], drop = True).where(formated.lat_f<=Config['filtering_parameters'][filter]['lat2'], drop = True)
            formated = formated.where(formated.depth_f>=Config['filtering_parameters'][filter]['bathy_min'], drop = True)
            formated = formated.where(formated.duration>=Config['filtering_parameters'][filter]['duration'], drop = True)
            formated.to_netcdf(Config['save_folder_path']+Config['tracking_formated_prefixe']+'Prop_'+Config['filtering_parameters'][filter]['name']+'_'+str(ind_depth)+'TEMP.nc')
            format = xr.open_dataset(Config['save_folder_path']+Config['tracking_formated_prefixe']+'Prop_'+Config['filtering_parameters'][filter]['name']+'_'+str(ind_depth)+'TEMP.nc')
            if len(format)>1 and len(format.eddy_num.where(format.time_i.dt.year>=Config['filtering_parameters'][filter]['yi'], drop = True).where(format.time_f.dt.year<=Config['filtering_parameters'][filter]['yf'], drop = True))>1:
                format = format.where(format.time_i.dt.year>=Config['filtering_parameters'][filter]['yi'], drop = True).where(format.time_f.dt.year<=Config['filtering_parameters'][filter]['yf'], drop = True)
                format['type'] = format['type'].astype(str)
            else:
                print('len tracked not good, depthi', ind_depth)
            format.to_netcdf(Config['save_folder_path']+Config['tracking_formated_prefixe']+'Prop_'+Config['filtering_parameters'][filter]['name']+'_'+str(ind_depth)+'_'+Config['suffixe_all']+'.nc')
            os.remove(Config['save_folder_path']+Config['tracking_formated_prefixe']+'Prop_'+Config['filtering_parameters'][filter]['name']+'_'+str(ind_depth)+'TEMP.nc')

            print('ind_depth ', ind_depth, ' done')
        else:
            print('ind_depth ', ind_depth, 'missing track')
            print(Config['save_folder_path'] + Config['prefixe_properties_tracking']+'tracked_depthi' + str(ind_depth) + Config['suffixe_all']+'tracked_properties.pickle')
    else:
        print('ind_depth ', ind_depth, ' already done')
    return None