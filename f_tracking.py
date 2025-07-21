import numpy as np
import xarray as xr
import os
from xorca.lib import load_xorca_dataset # type: ignore
import eddytools as et # type: ignore
from glob import glob
import pickle
import time
import cftime
import xgcm # type: ignore
import datetime
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.simplefilter('ignore', UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
import f_generic

def tracking_REF12(depth_index,yi,  yf, Config):
    domain_path = Config['domain_path']
    prefixe_interpolation = Config['prefixe_interpolation']
    prefixe_detection_m = Config['prefixe_detection_m']
    prefixe_tracking = Config['prefixe_tracking']
    save_folder_path = Config['save_folder_path']
    datapath = Config['datapath']
    meshpath = Config['meshpath']
    maskpath = Config['maskpath']
    tracking_parameters = Config['tracking_parameters']
    ds_zgr = xr.open_dataset(domain_path)
    L_depth_value = ds_zgr.nav_lev.values.astype(int)
    depth = L_depth_value[depth_index]
    tracking_parameters['data_path']= save_folder_path, # path to the detected eddies pickle files
    tracking_parameters['file_root']= prefixe_detection_m[:-1], 
    tracking_parameters['file_spec']= 'depthi' + str(depth_index) + Config['suffixe_all']+'_eddies',
    tracking_parameters['end_time']= str(yf)+'-12-31', # time range end
    tracking_parameters['start_time']= str(yi)+'-01-01', # time range start
    tracking_parameters['data_path'] = tracking_parameters['data_path'][0]
    tracking_parameters['file_root'] = tracking_parameters['file_root'][0]
    tracking_parameters['file_spec'] = tracking_parameters['file_spec'][0]
    tracking_parameters['end_time'] = tracking_parameters['end_time'][0]
    tracking_parameters['start_time'] = tracking_parameters['start_time'][0]
    # try tracking only if detection exists

    if not np.all([os.path.exists(Config['save_folder_path'] + Config['prefixe_detection_m'] + str(y)+'_'+'{0:02}'.format(m)+ '_depthi' + str(depth_index) + Config['suffixe_all']+'_eddies.pickle') for y in range(yi, yf+1) for m in range(1, 13)]): 
        print('I cannot find detection for that depth on all months ', str(yf)+'/12' )
    else:
        print('Tracking file', Config['save_folder_path'] + Config['prefixe_tracking']+str(depth_index)+Config['suffixe_all']+'_tracking.pickle')
        if not os.path.exists(Config['save_folder_path'] + Config['prefixe_tracking']+str(depth_index)+Config['suffixe_all']+'_tracking.pickle'):
            # Now we track the eddies, all information needed has to be added to `tracking_parameters`
            print(tracking_parameters)
            tracking = et.tracking.track(tracking_parameters, in_file=True)
            
            print('I found ', len(tracking), ' eddies')
            with open(Config['save_folder_path'] + Config['prefixe_tracking']+str(depth_index)+Config['suffixe_all']+'_tracking.pickle', 'wb') as f:
                pickle.dump(tracking, f, pickle.HIGHEST_PROTOCOL)
            f.close()        
        else:
            print('I already have the tracking for that depth ')
        get_depth_tracked(Config, yi, depth_index)    
    return None

def get_eddies_1d(depth_index, yi, yf, Config):
    file=Config['save_folder_path'] +Config['prefixe_tracking']+str(depth_index)+Config['suffixe_all']+'_tracking.pickle'
    infile = open(file,'rb')
    tracked = pickle.load(infile)
    infile.close()

    detected = []
    files=[Config['save_folder_path'] +Config['prefixe_detection_m']+str(y)+'_'+'{0:02}'.format(m)+'_depthi'+str(depth_index)+Config['suffixe_all']+'_eddies.pickle' for y in range(yi, yf+1) for m in range(1, 13)]
    c = 0
    if not np.all([os.path.exists(Config['save_folder_path'] + Config['prefixe_detection_m'] + str(y)+'_'+'{0:02}'.format(m)+ '_depthi' + str(depth_index) + Config['suffixe_all']+'_eddies.pickle') for y in range(yi, yf+1) for m in range(1, 13)]): 
        print('I cannot find detection for that depth on all months ', str(yf)+'/12' )
    else:  
        if not os.path.exists(Config['save_folder_path'] + Config['prefixe_tracking']+str(depth_index)+Config['suffixe_all']+'_tracking_1d.pickle'):
            for file in files:
                infile = open(file,'rb')
                detected.append(pickle.load(infile))
                infile.close()
                
            # get ID of tracked eddies
            In_Track = {}
            for etr in range(len(tracked)): #for each eddy tracked
                for tst in range(len(tracked[etr]['time'])): #for each time step of the tracked eddy
                    time_step = tracked[etr]['time'][tst]
                    day = (time_step- cftime.datetime(yi,1,1,12, calendar = 'noleap')).days
                    try:
                        In_Track[day] 
                    except:
                        In_Track[day] = {}
                    eddy_num = tracked[etr]['eddy_num'][tst]
                    try : 
                        In_Track[day][eddy_num]
                    except:
                        In_Track[day][eddy_num] = {}
                    In_Track[day][eddy_num]['InTrack'] = True
            # get characteristics of un-tracked eddies
            date_init = cftime.datetime(yi, 1, 1, 12, calendar = 'noleap')
            Tracked_dt_1 = {}
            c = 0
            for mo in range(len(files)):
                ddet = detected[mo]
                for da in range(len(ddet)):
                    det = ddet[da]
                    for e in range(len(det)):
                        #try:
                        i  = (det[e]['time'] - date_init).days
                        try:
                            In_Track[i][e]
                        except:
                            #print('Eddy ', e, 'on day ', d, 'was not tracked')
                            Tracked_dt_1[c] = {}
                            Tracked_dt_1[c]['time'] = det[e]['time']
                            Tracked_dt_1[c]['lon'] = det[e]['lon']
                            Tracked_dt_1[c]['lat'] = det[e]['lat']
                            Tracked_dt_1[c]['amp'] = det[e]['amp']
                            Tracked_dt_1[c]['eddy_j'] = det[e]['eddy_j']
                            Tracked_dt_1[c]['eddy_i'] = det[e]['eddy_i']
                            Tracked_dt_1[c]['area'] = det[e]['area']
                            Tracked_dt_1[c]['scale'] = det[e]['scale']
                            Tracked_dt_1[c]['type'] = det[e]['type']
                            Tracked_dt_1[c]['eddy_num'] = e
                            c+=1
                    #except:
                    #    print('I pass')
            # save 1d eddies properties
            with open(Config['save_folder_path'] + Config['prefixe_tracking']+str(depth_index)+Config['suffixe_all']+'_tracking_1d.pickle', 'wb') as f:
                pickle.dump(Tracked_dt_1, f, pickle.HIGHEST_PROTOCOL)
            f.close()    
        else:
            print('I already have tracking 1d at that depth')
        get_depth_tracked_1d(Config, yi, depth_index)
    return None

def get_depth_tracked(Config,yi, depth_index):
    file=Config['save_folder_path'] +Config['prefixe_tracking']+str(depth_index)+Config['suffixe_all']+'_tracking.pickle'
    infile = open(file,'rb')
    tracked = pickle.load(infile)
    infile.close()
    try:
        tracked[0]['depth']
        print('I already have the depth tracked at that depth !')
    except:
        print('145 : domain path', Config['domain_path'])
        ds_zgr = xr.open_dataset(Config['domain_path'])
        print('148 : xorca', Config['one_T_file_path'])
        dsT = load_xorca_dataset(data_files=[Config['one_T_file_path']],model_config='NEST', aux_files=[Config['domain_path']])
        dsT['at'] = dsT['e1t']*dsT['e2t']
        dsT['au'] = dsT['e1u']*dsT['e2u']
        dsT['av'] = dsT['e1v']*dsT['e2v']
        dsT['af'] = dsT['e1f']*dsT['e2f']
        dsT = dsT.set_coords(('at', 'au', 'av', 'af')).isel(z_c = 0)
        dsT['bathy_meter'] = (('y_c', 'x_c'), ds_zgr.isel(t = 0).bathy_meter.values)
        metrics = {('X',): ['e1u', 'e1v'], ('Y',): ['e2v', 'e2u'], ('X', 'Y'): ['at', 'au', 'av', 'af']}
        coords={"X": {"center": "x_c", "right": "x_r"}, "Y":{'center':'y_c', 'right':'y_r'}}
        grid = xgcm.Grid(dsT, periodic=["X"], coords = coords, metrics=metrics, boundary = 'fill', fill_value=np.nan)
        print('159 : savegrid')
        ds_grid = f_generic.save_grid_REF12(Config)

        def get_ind_env(eddy_i_ind, eddy_j_ind, fact, j_max = dsT.y_r.max().values, i_max = dsT.x_r.max().values):
            eddy_i_size = (eddy_i_ind.max() - eddy_i_ind.min()).values
            eddy_j_size = (eddy_j_ind.max() - eddy_j_ind.min()).values
            eddy_j_aug = np.arange(eddy_j_ind.min().values - fact*eddy_j_size, eddy_j_ind.max().values + fact*eddy_j_size +1, 1)
            eddy_i_aug = np.arange(eddy_i_ind.min().values - fact*eddy_i_size, eddy_i_ind.max().values + fact*eddy_i_size +1, 1)

            II, JJ = np.meshgrid(eddy_i_aug, eddy_j_aug)
            II = II.flatten(); JJ = JJ.flatten()
            L = [(x,y) for (x,y) in zip(II,JJ)]

            for (i,j) in zip(eddy_i_ind.values,eddy_j_ind.values):
                L.remove((i,j))
            L = np.array(L)
            J = L[:,1]
            I = L[:,0]
            J[np.where(J>j_max)[0]] = 2*j_max - J[np.where(J>j_max)[0]]
            I[np.where(J>j_max)[0]] = i_max - I[np.where(J>j_max)[0]]
            return I,J
        print('180 : bathy interp')
        Bathy_f = grid.interp(dsT.bathy_meter, axis = ('X', 'Y')).compute()

        for tr in range(len(tracked)):
            if tr%int(len(tracked)/100)==0:
                print(np.round(tr/len(tracked)*100, 1), '%')
            e = tracked[tr]
            depth_t = np.zeros((len(e['time'])))
            for t in range(len(e['time'])):
                eddy_j_ind = ds_grid.y[e['eddy_j'][t]]
                eddy_i_ind = ds_grid.x[e['eddy_i'][t]]
                depth_t[t] = Bathy_f.sel(x_r = eddy_i_ind).sel(y_r = eddy_j_ind).mean()
            e['depth'] = depth_t
        print('194 : save new tracked')
        # save into tracked 
        with open(Config['save_folder_path'] + Config['prefixe_tracking']+str(depth_index)+Config['suffixe_all']+'_tracking.pickle', 'wb') as f:
            pickle.dump(tracked, f, pickle.HIGHEST_PROTOCOL)
        f.close()    
    return None

def get_depth_tracked_1d(Config, yi, depth_index):
    file=Config['save_folder_path'] +Config['prefixe_tracking']+str(depth_index)+Config['suffixe_all']+'_tracking_1d.pickle'
    infile = open(file,'rb')
    tracked = pickle.load(infile)
    infile.close()
    try:
        tracked[0]['depth']
        print('I already have the depth tracked at that depth !')
    except:
        ds_zgr = xr.open_dataset(Config['domain_path'])
        dsT = load_xorca_dataset(data_files=[Config['one_T_file_path']],model_config='NEST', aux_files=[Config['domain_path']])
        dsT['at'] = dsT['e1t']*dsT['e2t']
        dsT['au'] = dsT['e1u']*dsT['e2u']
        dsT['av'] = dsT['e1v']*dsT['e2v']
        dsT['af'] = dsT['e1f']*dsT['e2f']
        dsT = dsT.set_coords(('at', 'au', 'av', 'af')).isel(z_c = 0)
        dsT['bathy_meter'] = (('y_c', 'x_c'), ds_zgr.isel(t = 0).bathy_meter.values)
        metrics = {('X',): ['e1u', 'e1v'], ('Y',): ['e2v', 'e2u'], ('X', 'Y'): ['at', 'au', 'av', 'af']}
        coords={"X": {"center": "x_c", "right": "x_r"}, "Y":{'center':'y_c', 'right':'y_r'}}
        grid = xgcm.Grid(dsT, periodic=["X"], coords = coords, metrics=metrics, boundary = 'fill', fill_value=np.nan)
        ds_grid = f_generic.save_grid_REF12(Config)

        def get_ind_env(eddy_i_ind, eddy_j_ind, fact, j_max = dsT.y_r.max().values, i_max = dsT.x_r.max().values):
            eddy_i_size = (eddy_i_ind.max() - eddy_i_ind.min()).values
            eddy_j_size = (eddy_j_ind.max() - eddy_j_ind.min()).values
            eddy_j_aug = np.arange(eddy_j_ind.min().values - fact*eddy_j_size, eddy_j_ind.max().values + fact*eddy_j_size +1, 1)
            eddy_i_aug = np.arange(eddy_i_ind.min().values - fact*eddy_i_size, eddy_i_ind.max().values + fact*eddy_i_size +1, 1)

            II, JJ = np.meshgrid(eddy_i_aug, eddy_j_aug)
            II = II.flatten(); JJ = JJ.flatten()
            L = [(x,y) for (x,y) in zip(II,JJ)]

            for (i,j) in zip(eddy_i_ind.values,eddy_j_ind.values):
                L.remove((i,j))
            L = np.array(L)
            J = L[:,1]
            I = L[:,0]
            J[np.where(J>j_max)[0]] = 2*j_max - J[np.where(J>j_max)[0]]
            I[np.where(J>j_max)[0]] = i_max - I[np.where(J>j_max)[0]]
            return I,J

        Bathy_f = grid.interp(dsT.bathy_meter, axis = ('X', 'Y')).compute()

        for tr in range(len(tracked)):
            if tr%int(len(tracked)/100)==0:
                print(np.round(tr/len(tracked)*100, 1), '%')

            e = tracked[tr]

            eddy_j_ind = ds_grid.y[e['eddy_j']]
            eddy_i_ind = ds_grid.x[e['eddy_i']]
            depth_t = Bathy_f.sel(x_r = eddy_i_ind).sel(y_r = eddy_j_ind).mean().values
            e['depth'] = depth_t

        # save into tracked
        with open(Config['save_folder_path'] + Config['prefixe_tracking']+str(depth_index)+Config['suffixe_all']+'_tracking_1d.pickle', 'wb') as f:
            pickle.dump(tracked, f, pickle.HIGHEST_PROTOCOL)
        f.close()    
    return None    

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

def track_properties_REF12(depth_ind, Config, yi, yf):
    if np.all(np.array([os.path.exists(Config['save_folder_path']+Config['prefixe_properties_m']+str(y)+'_{0:02}'.format(m)+ '_depthi' + str(depth_ind) + Config['suffixe_all']+'.pickle') for y in range(yi, yf+1) for m in range(1, 13)])) and np.all(np.array([os.path.getsize(Config['save_folder_path']+Config['prefixe_properties_m']+str(y)+'_{0:02}'.format(m)+ '_depthi' + str(depth_ind) + Config['suffixe_all']+'.pickle')>10000 for y in range(yi, yf+1) for m in range(1, 13)])):

        Properties_extract = Config['properties_parameters']['Properties_extract']
        #---
        tr = load_tracked(Config['save_folder_path'], depth_ind ,Config['prefixe_tracking'], Config['suffixe_all'])
        print(Config['save_folder_path'] + Config['prefixe_properties_tracking']+'tracked_depthi' + str(depth_ind) + Config['suffixe_all']+'tracked_properties.pickle')
        print('load pr in process')
        if os.path.exists(Config['save_folder_path'] + Config['prefixe_properties_tracking']+'tracked_depthi' + str(depth_ind) + Config['suffixe_all']+'tracked_properties.pickle'):       
            infile = open(Config['save_folder_path'] + Config['prefixe_properties_tracking']+'tracked_depthi' + str(depth_ind) + Config['suffixe_all']+'tracked_properties.pickle','rb')
            pr = pickle.load(infile)
            infile.close()
        else:
            print('reset pr')
            pr = {}
        
        for te in range(len(pr), len(tr)):
            if te%int(len(tr)/1000)==0:
                print(np.round(te/len(tr)*100, 1), '%')
            timeframes = tr[te]['time']
            y0 = timeframes[0].year
            m0 = timeframes[0].month
            infile = open(Config['save_folder_path'] + Config['prefixe_properties_m'] + str(y0)+'_'+'{0:02}'.format(m0) + '_depthi' + str(depth_ind) + Config['suffixe_all']+'.pickle','rb' )
            properties = pickle.load(infile)
            infile.close()    
            eddy_nums = tr[te]['eddy_num']
            Pr_t = {}
            Properties = np.nan*np.zeros((len(Properties_extract), len(timeframes)))

            for t in range(len(timeframes)):
                if timeframes[t].year != y0 or timeframes[t].month != m0:
                    y0 = timeframes[t].year
                    m0 = timeframes[t].month
                    try:
                        with open(Config['save_folder_path'] + Config['prefixe_properties_m'] + str(y0)+'_'+'{0:02}'.format(m0) + '_depthi' + str(depth_ind) + Config['suffixe_all']+'.pickle','rb' ) as infile:
                            properties = pickle.load(infile)
                        infile.close()
                    except:
                        print(Config['save_folder_path'] + Config['prefixe_properties_m'] + str(y0)+'_'+'{0:02}'.format(m0) + '_depthi' + str(depth_ind) + Config['suffixe_all']+'.pickle','rb' )

                for ip in range(len(Properties_extract)):
                    Properties[ip,t] = float(properties[timeframes[t].day-1][eddy_nums[t]][Properties_extract[ip]])

            for ip in range(len(Properties_extract)):
                p = Properties_extract[ip]
                Pr_t[p] = Properties[ip,:]
            Pr_t['time'] = tr[te]['time']
            Pr_t['lon'] = tr[te]['lon']
            Pr_t['lat'] = tr[te]['lat']
            del Properties
            pr[te] = Pr_t
            if te%int(len(tr)/100)==0:
                print('Saving', np.round(te/len(tr)*100, 1), '%')
                with open(Config['save_folder_path'] + Config['prefixe_properties_tracking']+'tracked_depthi' + str(depth_ind) + Config['suffixe_all']+'tracked_properties.pickle', 'wb') as f:
                    pickle.dump(pr, f, pickle.HIGHEST_PROTOCOL)
                f.close()
        with open(Config['save_folder_path'] + Config['prefixe_properties_tracking']+'tracked_depthi' + str(depth_ind) + Config['suffixe_all']+'tracked_properties.pickle', 'wb') as f:
            pickle.dump(pr, f, pickle.HIGHEST_PROTOCOL)
        f.close()
    else:
        print('I miss some properties')
        for y in range(yi, yf+1):
            for m in range(1,13):
                if os.path.exists(Config['save_folder_path']+Config['prefixe_properties_m']+str(y)+'_{0:02}'.format(m)+ '_depthi' + str(depth_ind) + Config['suffixe_all']+'.pickle'):
                    if os.path.getsize(Config['save_folder_path']+Config['prefixe_properties_m']+str(y)+'_{0:02}'.format(m)+ '_depthi' + str(depth_ind) + Config['suffixe_all']+'.pickle')<10000:
                        print('file is too small',y, '/', m, ' : ', os.path.getsize(Config['save_folder_path']+Config['prefixe_properties_m']+str(y)+'_{0:02}'.format(m)+ '_depthi' + str(depth_ind) + Config['suffixe_all']+'.pickle'))
                else:
                    print('I miss : ', y, '/', m)
    return None    

def track_properties_REF12_1d(depth_ind, Config, yi, yf):
    if np.all(np.array([os.path.exists(Config['save_folder_path']+Config['prefixe_properties_m']+str(y)+'_{0:02}'.format(m)+ '_depthi' + str(depth_ind) + Config['suffixe_all']+'.pickle') for y in range(yi, yf+1) for m in range(1, 13)])) and np.all(np.array([os.path.getsize(Config['save_folder_path']+Config['prefixe_properties_m']+str(y)+'_{0:02}'.format(m)+ '_depthi' + str(depth_ind) + Config['suffixe_all']+'.pickle')>10000 for y in range(yi, yf+1) for m in range(1, 13)])):

        Properties_extract = Config['properties_parameters']['Properties_extract']
        #---
        tr = load_tracked_1d(Config['save_folder_path'], depth_ind ,Config['prefixe_tracking'], Config['suffixe_all'])
        if os.path.exists(Config['save_folder_path'] + Config['prefixe_properties_tracking']+'tracked_depthi' + str(depth_ind) + Config['suffixe_all']+'tracked_properties_1d.pickle'):       
            infile = open(Config['save_folder_path'] + Config['prefixe_properties_tracking']+'tracked_depthi' + str(depth_ind) + Config['suffixe_all']+'tracked_properties_1d.pickle','rb')
            pr = pickle.load(infile)
            infile.close()
        else:
            print('reset pr')
            pr = {}
        
        for te in range(len(pr), len(tr)):
            if te%int(len(tr)/1000)==0:
                print(np.round(te/len(tr)*100, 1), '%')
            timeframes = tr[te]['time']
            infile = open(Config['save_folder_path'] + Config['prefixe_properties_m'] + str(timeframes)[:4]+'_'+str(timeframes)[5:7] + '_depthi' + str(depth_ind) + Config['suffixe_all']+'.pickle','rb' )
            properties = pickle.load(infile)
            infile.close()    
            eddy_nums = tr[te]['eddy_num']
            Pr_t = {}
            Properties = np.nan*np.zeros((len(Properties_extract)))
            for ip in range(len(Properties_extract)):
                p = Properties_extract[ip]
                Pr_t[p] = float(properties[int(str(timeframes)[8:10])-1][eddy_nums][Properties_extract[ip]])
            Pr_t['time'] = tr[te]['time']
            Pr_t['lon'] = tr[te]['lon']
            Pr_t['lat'] = tr[te]['lat']
            del Properties
            pr[te] = Pr_t


            if te%int(len(tr)/100)==0:
                print('Saving', np.round(te/len(tr)*100, 1), '%')
                with open(Config['save_folder_path'] + Config['prefixe_properties_tracking']+'tracked_depthi' + str(depth_ind) + Config['suffixe_all']+'tracked_properties_1d.pickle', 'wb') as f:
                    pickle.dump(pr, f, pickle.HIGHEST_PROTOCOL)
                f.close()
        with open(Config['save_folder_path'] + Config['prefixe_properties_tracking']+'tracked_depthi' + str(depth_ind) + Config['suffixe_all']+'tracked_properties_1d.pickle', 'wb') as f:
            pickle.dump(pr, f, pickle.HIGHEST_PROTOCOL)
        f.close()
    else:
        print('I miss some properties')
        for y in range(yi, yf+1):
            for m in range(1,13):
                if os.path.exists(Config['save_folder_path']+Config['prefixe_properties_m']+str(y)+'_{0:02}'.format(m)+ '_depthi' + str(depth_ind) + Config['suffixe_all']+'.pickle'):
                    if os.path.getsize(Config['save_folder_path']+Config['prefixe_properties_m']+str(y)+'_{0:02}'.format(m)+ '_depthi' + str(depth_ind) + Config['suffixe_all']+'.pickle')<10000:
                        print('file is too small', y, '/', m, ' : ', os.path.getsize(Config['save_folder_path']+Config['prefixe_properties_m']+str(y)+'_{0:02}'.format(m)+ '_depthi' + str(depth_ind) + Config['suffixe_all']+'.pickle'))
                else:
                    print('I miss : ', y, '/', m)
    return None    

