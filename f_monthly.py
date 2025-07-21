from xorca.lib import load_xorca_dataset # type: ignore
import numpy as np
import xarray as xr
import os
import pickle
import cftime
import gsw
import datetime
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.simplefilter('ignore', UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
import f_generic


air_frac = 0

def make_ORCA_monthly_m_filter_dt(yi, yf, depth_ind, Config, filter_time = 2):
    save_folder_path = Config['save_folder_path']
    domain_path = Config['domain_path']
    ds_zgr = xr.open_dataset(domain_path)
    L_depth_value = ds_zgr.nav_lev.values.astype(int)
    depth = L_depth_value[depth_ind]

    # empty array
    de = f_generic.do_empty_array_REF12_fpoints(Config)

    ds_grid = f_generic.save_grid_REF12(Config)

    COORDx = np.zeros((de.x.shape[0], ))
    for i in range(ds_grid.x.shape[0]):
        ind = np.where(de.x.values +0.5 == ds_grid.x.values[i])[0]
        COORDx[i] = ind
        
    COORDy = np.zeros((de.y.shape[0], ))
    for i in range(ds_grid.y.shape[0]):
        ind = np.where(de.y.values+0.5 == ds_grid.y.values[i])[0]
        COORDy[i] = ind
    infile = open(Config['save_folder_path'] + Config['prefixe_tracking']  +str(depth_ind) + Config['suffixe_all']+'_tracking.pickle','rb')
    tracked = pickle.load(infile)
    infile.close()
    #loop month
    for y in range(yi, yf+1):
        for m in range(1, 13):
            print(y, '/', m, '-', depth)
            if not os.path.exists(Config['save_folder_path']+Config['prefixe_monthy_means']+str(y)+'_'+'{0:02}'.format(m)+'_depthi'+str(depth_ind)+'_'+Config['suffixe_all']+'.nc'):
                get_tracking_info_REF12_m(yi, yf, depth_ind, Config)

                dsAll = de.copy(deep = True)
                date = datetime.datetime(year = y, month = m, day = 15)
                dsAll['time_counter'] = (('t'), np.array([date]))
                dsAll = dsAll.set_coords('time_counter')
                dsAll = dsAll.swap_dims({'t':'time_counter'})
                dsAll = dsAll.isel(z = depth_ind)
                try:
                    dsAll = dsAll.drop({'t'})
                except:
                    pass

                try:
                    nmb_days = (cftime.datetime(y,m,31, calendar = 'NoLeap')-cftime.datetime(y,m,1, calendar = 'NoLeap')).days+1
                except:
                    try:
                        nmb_days = (cftime.datetime(y,m,30, calendar = 'NoLeap')-cftime.datetime(y,m,1, calendar = 'NoLeap')).days+1
                    except:
                        nmb_days = (cftime.datetime(y,m,28, calendar = 'NoLeap')-cftime.datetime(y,m,1, calendar = 'NoLeap')).days+1

                Nb = np.zeros((nmb_days, de.y.shape[0], de.x.shape[0]))
                NCb = np.zeros((nmb_days, de.y.shape[0], de.x.shape[0]))
                Rb = np.nan*np.zeros((nmb_days, de.y.shape[0], de.x.shape[0]))
                Ab = np.nan*np.zeros((nmb_days, de.y.shape[0], de.x.shape[0]))
                dtb = np.nan*np.zeros((nmb_days, de.y.shape[0], de.x.shape[0]))
                D = np.nan*np.zeros((nmb_days, de.y.shape[0], de.x.shape[0]))
                time_step = np.nan*np.zeros((nmb_days, de.y.shape[0], de.x.shape[0])) #birth at time step 0

                infile = open(Config['save_folder_path']+Config['prefixe_detection_m']+str(y)+'_'+'{0:02}'.format(m)+'_depthi'+str(depth_ind)+Config['suffixe_all']+'_eddies.pickle','rb')
                DET = pickle.load(infile)
                infile.close()
                infile = open(Config['save_folder_path'] + Config['prefixe_tracking_det_indexes_m'] + str(y)+'-'+'{0:02}'.format(m) + '_depthi' + str(depth_ind) + Config['suffixe_all']+'_track_det.pickle','rb')
                TRACK = pickle.load(infile)
                infile.close()  
                for t in range(nmb_days):
                    print(t, '/', nmb_days)
                    timeframe = cftime.datetime(y, m, 1, calendar = 'NoLeap') + datetime.timedelta(days = t)
                    data = DET[t]
                    tracking_info = TRACK[t]
                    detection_exist = len(data)>0
                    if detection_exist:
                        if str(data[0]['time'])[:10] != str(timeframe)[0:10]:
                            print('I have a time problem with depth ', depth, 'at time ', timeframe)
                            break
                        else:
                            for tr in tracking_info:
                                time_step[t,COORDy[data[tr]['eddy_j']].astype(int),COORDx[data[tr]['eddy_i']].astype(int)] = tracking_info[tr]['time_step']
                                if  tracking_info[tr]['time_step']==0:
                                    # birth values only
                                    Nb[t,COORDy[data[tr]['eddy_j']].astype(int),COORDx[data[tr]['eddy_i']].astype(int)] = 1
                                    Rb[t,COORDy[data[tr]['eddy_j']].astype(int),COORDx[data[tr]['eddy_i']].astype(int)] = data[tr]['scale']
                                    Ab[t,COORDy[data[tr]['eddy_j']].astype(int),COORDx[data[tr]['eddy_i']].astype(int)] = data[tr]['amp']
                                    if data[tr]['type'] =='cyclonic':
                                        NCb[t,COORDy[data[tr]['eddy_j']].astype(int),COORDx[data[tr]['eddy_i']].astype(int)] = 1

                                    dtb[t,COORDy[data[tr]['eddy_j']].astype(int),COORDx[data[tr]['eddy_i']].astype(int)] = tracking_info[tr]['dt']
                                if  tracking_info[tr]['time_step']==0 and tracking_info[tr]['track_coord']<len(tracked):
                                    # birth values only
                                    D[t,COORDy[data[tr]['eddy_j']].astype(int),COORDx[data[tr]['eddy_i']].astype(int)] = f_generic.get_distance_latlon2m(tracked[tracking_info[tr]['track_coord']]['lon'][0], tracked[tracking_info[tr]['track_coord']]['lat'][0], tracked[tracking_info[tr]['track_coord']]['lon'][-1], tracked[tracking_info[tr]['track_coord']]['lat'][-1])
                                
                    else:
                        print('No detection on day ', timeframe)

                dstemp = de.copy(deep = True)
                dstemp = dstemp.squeeze()
                dstemp['day'] = (('day'), np.arange(nmb_days))
                dstemp = dstemp.set_coords('day')
                dstemp.drop('time_counter')

                dstemp['Nb'] = (('day', 'y', 'x'), Nb)
                dstemp['NCb'] = (('day', 'y', 'x'), NCb)
                dstemp['Rb'] = (('day', 'y', 'x'), Rb)
                dstemp['Ab'] = (('day', 'y', 'x'), Ab)
                dstemp['dtb'] = (('day', 'y', 'x'), dtb)
                dstemp['D'] = (('day', 'y', 'x'), D)
                dstemp['timestep'] = (('day', 'y', 'x'), time_step)



                dsAll['Nb_sum'] = dstemp.Nb.where(dstemp.dtb>=filter_time).sum('day', skipna = True).astype(np.uint16).expand_dims('time_counter')

                dsAll['NCb_sum'] = dstemp.NCb.where(dstemp.dtb>=filter_time).sum('day').expand_dims( 'time_counter')
                
                dsAll['Rb_mean'] = dstemp.Rb.where(dstemp.dtb>=filter_time).mean('day', skipna = True).expand_dims( 'time_counter')

                dsAll['Ab_mean'] = dstemp.Ab.where(dstemp.dtb>=filter_time).mean('day', skipna = True).expand_dims( 'time_counter')

                dsAll['dtb_mean'] = dstemp.dtb.where(dstemp.dtb>=filter_time).mean('day', skipna = True).expand_dims( 'time_counter')

                dsAll['time_step_mean'] = dstemp.timestep.where(dstemp.dtb>=filter_time).mean('day', skipna = True).expand_dims( 'time_counter')
                
                dsAll['D'] = dstemp.D.where(dstemp.dtb>=filter_time).mean('day', skipna = True).expand_dims( 'time_counter')

                #save ds0 y-m
                dsAll.to_netcdf(Config['save_folder_path']+Config['prefixe_monthy_means']+str(y)+'_'+'{0:02}'.format(m)+'_depthi'+str(depth_ind)+'_'+Config['suffixe_all']+'.nc')
            else:
                print('I already have dsAll for y: ', y, ' m : ', m)
    return None

def get_tracking_info_REF12_m(yi, yf, depth_ind, Config):       
    save_folder_path = Config['save_folder_path']
    prefixe_tracking = Config['prefixe_tracking']
    suffixe_all = Config['suffixe_all']
    Det_Index = {}
    date_init = cftime.datetime(yi, 1, 1, calendar=u'NoLeap')
    date_fin = cftime.datetime(yf, 12, 31, calendar=u'NoLeap')
    I_have_tracks = os.path.exists(save_folder_path+prefixe_tracking+str(depth_ind)+suffixe_all+'_tracking.pickle') and os.path.exists(save_folder_path+prefixe_tracking+str(depth_ind)+suffixe_all+'_tracking_1d.pickle')
    if not I_have_tracks:
        print('I need both trackings 1d and 2d for this step')
    else:
        if not os.path.exists(save_folder_path + Config['prefixe_tracking_det_indexes_m'] + str(yf)+'-'+'{0:02}'.format(12)+ '_depthi' + str(depth_ind) + suffixe_all+'_track_det.pickle'):
            infile = open(save_folder_path+prefixe_tracking+str(depth_ind)+suffixe_all+'_tracking.pickle','rb')
            tracking = pickle.load(infile)
            infile.close()
            infile = open(save_folder_path+prefixe_tracking+str(depth_ind)+suffixe_all+'_tracking_1d.pickle','rb')
            tracking1d = pickle.load(infile)
            infile.close()
            for ind_tr in range(len(tracking)):
                ed_tracked = tracking[ind_tr]
                for ind_time in range(len(ed_tracked['time'])):
                    time = ed_tracked['time'][ind_time]
                    eddy_num = ed_tracked['eddy_num'][ind_time]
                    i = (time - date_init).days
                    try:
                        Det_Index[i] 
                    except:
                        Det_Index[i] = {}
                    try: 
                        Det_Index[i][eddy_num]
                    except:
                        Det_Index[i][eddy_num] = {}
                    Det_Index[i][eddy_num]['track_coord'] = ind_tr
                    Det_Index[i][eddy_num]['dt'] = len(ed_tracked['time'])
                    Det_Index[i][eddy_num]['time_step'] = ind_time
                    
            for ind_tr in range(len(tracking1d)):
                ed_tracked = tracking1d[ind_tr]
                time = ed_tracked['time']
                eddy_num = ed_tracked['eddy_num']
                i = (time - date_init).days
                try:
                    Det_Index[i] 
                except:
                    Det_Index[i] = {}
                try: 
                    Det_Index[i][eddy_num]
                    print('This eddy was already detected - why ?', i, eddy_num, time, ind_tr)
                except:
                    Det_Index[i][eddy_num] = {}
                    pass
                Det_Index[i][eddy_num]['track_coord'] = ind_tr
                Det_Index[i][eddy_num]['dt'] = 1
                Det_Index[i][eddy_num]['time_step'] = 0
                i=0
            for y in range(yi, yf+1):
                for m in range(1, 13):
                    try:
                        df = cftime.datetime(y, m, 31, calendar=u'NoLeap')
                    except:
                        try:
                            df = cftime.datetime(y, m, 30, calendar=u'NoLeap')
                        except:
                            df = cftime.datetime(y, m, 28, calendar=u'NoLeap')
                    di = cftime.datetime(y, m, 1, calendar=u'NoLeap')
                    Det_track_m = {}
                    for j in range(((df-di)).days+1):
                        timeframe = di +datetime.timedelta(days = j)
                        try:
                            Det_track_m[j] = Det_Index[i]
                            i+=1
                        except:
                            Det_track_m[j] = {}
                            i+=1
                            print('I miss ', timeframe, ', no eddies tracked on that day ')
                    with open(save_folder_path + Config['prefixe_tracking_det_indexes_m'] + str(y)+'-'+'{0:02}'.format(m)+ '_depthi' + str(depth_ind) + suffixe_all+'_track_det.pickle', 'wb') as f:
                        pickle.dump(Det_track_m, f, pickle.HIGHEST_PROTOCOL)
                    f.close()
    return None

def make_ORCA_monthly_m_filter_prop(yi, yf, depth_ind, Config, filter_time = 2):
    save_folder_path = Config['save_folder_path']
    domain_path = Config['domain_path']
    ds_zgr = xr.open_dataset(domain_path)
    L_depth_value = ds_zgr.nav_lev.values.astype(int)
    depth = L_depth_value[depth_ind]

    # empty array
    de = f_generic.do_empty_array_REF12_fpoints(Config)

    ds_grid = f_generic.save_grid_REF12(Config)

    COORDx = np.zeros((de.x.shape[0], ))
    for i in range(ds_grid.x.shape[0]):
        ind = np.where(de.x.values +0.5 == ds_grid.x.values[i])[0]
        COORDx[i] = ind
        
    COORDy = np.zeros((de.y.shape[0], ))
    for i in range(ds_grid.y.shape[0]):
        ind = np.where(de.y.values+0.5 == ds_grid.y.values[i])[0]
        COORDy[i] = ind
    infile = open(Config['save_folder_path'] + Config['prefixe_tracking']  +str(depth_ind) + Config['suffixe_all']+'_tracking.pickle','rb')
    tracked = pickle.load(infile)
    infile.close()
    #loop month
    for y in range(yi, yf+1):
        for m in range(1, 13):
            print(y, '/', m, '-', depth)
            if not os.path.exists(Config['save_folder_path']+Config['prefixe_monthyP_means']+str(y)+'_'+'{0:02}'.format(m)+'_depthi'+str(depth_ind)+'_'+Config['suffixe_all']+'.nc'):
                get_tracking_info_REF12_m(yi, yf, depth_ind, Config)

                dsAll = de.copy(deep = True)
                date = datetime.datetime(year = y, month = m, day = 15)
                dsAll['time_counter'] = (('t'), np.array([date]))
                dsAll = dsAll.set_coords('time_counter')
                dsAll = dsAll.swap_dims({'t':'time_counter'})
                dsAll = dsAll.isel(z = depth_ind)
                try:
                    dsAll = dsAll.drop({'t'})
                except:
                    pass

                try:
                    nmb_days = (cftime.datetime(y,m,31, calendar = 'NoLeap')-cftime.datetime(y,m,1, calendar = 'NoLeap')).days+1
                except:
                    try:
                        nmb_days = (cftime.datetime(y,m,30, calendar = 'NoLeap')-cftime.datetime(y,m,1, calendar = 'NoLeap')).days+1
                    except:
                        nmb_days = (cftime.datetime(y,m,28, calendar = 'NoLeap')-cftime.datetime(y,m,1, calendar = 'NoLeap')).days+1
        
                
                Nb = np.zeros((nmb_days, de.y.shape[0], de.x.shape[0]))
                T_m_i = np.nan*np.zeros((nmb_days, de.y.shape[0], de.x.shape[0]))
                S_m_i = np.nan*np.zeros((nmb_days, de.y.shape[0], de.x.shape[0]))
                DT_me_i = np.nan*np.zeros((nmb_days, de.y.shape[0], de.x.shape[0]))
                DS_me_i = np.nan*np.zeros((nmb_days, de.y.shape[0], de.x.shape[0]))
                dT_me_i = np.nan*np.zeros((nmb_days, de.y.shape[0], de.x.shape[0]))
                dS_me_i = np.nan*np.zeros((nmb_days, de.y.shape[0], de.x.shape[0]))
                dTf_me_i = np.nan*np.zeros((nmb_days, de.y.shape[0], de.x.shape[0]))
                dtb = np.nan*np.zeros((nmb_days, de.y.shape[0], de.x.shape[0]))

                with open(Config['save_folder_path']+Config['prefixe_detection_m']+str(y)+'_'+'{0:02}'.format(m)+'_depthi'+str(depth_ind)+Config['suffixe_all']+'_eddies.pickle','rb') as infile:
                    DET = pickle.load(infile)
                infile.close()
                with open(Config['save_folder_path']+Config['prefixe_properties_m']+str(y)+'_'+'{0:02}'.format(m)+ '_depthi' + str(depth_ind) + Config['suffixe_all']+'.pickle','rb') as infile:
                    PROP = pickle.load(infile)
                infile.close()
                with open(Config['save_folder_path'] + Config['prefixe_tracking_det_indexes_m'] + str(y)+'-'+'{0:02}'.format(m) + '_depthi' + str(depth_ind) + Config['suffixe_all']+'_track_det.pickle','rb') as infile:
                    TRACK = pickle.load(infile)
                infile.close()  
                for t in range(nmb_days):
                    print(t, '/', nmb_days)
                    timeframe = cftime.datetime(y, m, 1, calendar = 'NoLeap') + datetime.timedelta(days = t)
                    data = DET[t]
                    prop = PROP[t]
                    tracking_info = TRACK[t]
                    detection_exist = len(data)>0
                    if detection_exist:
                        if str(data[0]['time'])[:10] != str(timeframe)[0:10]:
                            print('I have a time problem with depth ', depth, 'at time ', timeframe)
                            break
                        else:
                            for tr in tracking_info:
                                if  tracking_info[tr]['time_step']==0:
                                    # birth values only
                                    Nb[t,COORDy[data[tr]['eddy_j']].astype(int),COORDx[data[tr]['eddy_i']].astype(int)] = 1
                                    T_m_i[t,COORDy[data[tr]['eddy_j']].astype(int),COORDx[data[tr]['eddy_i']].astype(int)] = prop[tr]['votemper_m']
                                    S_m_i[t,COORDy[data[tr]['eddy_j']].astype(int),COORDx[data[tr]['eddy_i']].astype(int)] = prop[tr]['vosaline_m']
                                    DT_me_i[t,COORDy[data[tr]['eddy_j']].astype(int),COORDx[data[tr]['eddy_i']].astype(int)] = prop[tr]['votemper_m'] - gsw.CT_freezing(prop[tr]['vosaline_m'],depth, air_frac) - prop[tr]['votemper_e'] +gsw.CT_freezing(prop[tr]['vosaline_e'],depth, air_frac)
                                    DS_me_i[t,COORDy[data[tr]['eddy_j']].astype(int),COORDx[data[tr]['eddy_i']].astype(int)] = prop[tr]['vosaline_m'] - prop[tr]['vosaline_e']
                                    dtb[t,COORDy[data[tr]['eddy_j']].astype(int),COORDx[data[tr]['eddy_i']].astype(int)] = tracking_info[tr]['dt']
                                    if prop[tr]['vosaline_m'] - prop[tr]['vosaline_e'] >0:
                                        dS_me_i[t,COORDy[data[tr]['eddy_j']].astype(int),COORDx[data[tr]['eddy_i']].astype(int)] = prop[tr]['vosaline_m'] - prop[tr]['vosaline_e'] - prop[tr]['vosaline_std_e']
                                    else:
                                        dS_me_i[t,COORDy[data[tr]['eddy_j']].astype(int),COORDx[data[tr]['eddy_i']].astype(int)] = prop[tr]['vosaline_m'] - prop[tr]['vosaline_e'] + prop[tr]['vosaline_std_e']
                                    if prop[tr]['votemper_m'] - gsw.CT_freezing(prop[tr]['vosaline_m'],depth, air_frac) - prop[tr]['votemper_e'] +gsw.CT_freezing(prop[tr]['vosaline_e'],depth, air_frac) >0:
                                        dT_me_i[t,COORDy[data[tr]['eddy_j']].astype(int),COORDx[data[tr]['eddy_i']].astype(int)] = prop[tr]['votemper_m'] - gsw.CT_freezing(prop[tr]['vosaline_m'],depth, air_frac) - prop[tr]['votemper_e'] +gsw.CT_freezing(prop[tr]['vosaline_e'],depth, air_frac) - prop[tr]['votemper_std_e']
                                        dTf_me_i[t,COORDy[data[tr]['eddy_j']].astype(int),COORDx[data[tr]['eddy_i']].astype(int)] = prop[tr]['votemper_m'] - gsw.CT_freezing(prop[tr]['vosaline_m'],depth, air_frac) - prop[tr]['votemper_e'] +gsw.CT_freezing(prop[tr]['vosaline_e'],depth, air_frac) - prop[tr]['votemper_std_e'] - gsw.CT_freezing(prop[tr]['vosaline_e']-prop[tr]['vosaline_std_e'],depth, air_frac)+gsw.CT_freezing(prop[tr]['vosaline_e']+prop[tr]['vosaline_std_e'],depth, air_frac)
                                    else :
                                        dT_me_i[t,COORDy[data[tr]['eddy_j']].astype(int),COORDx[data[tr]['eddy_i']].astype(int)] = prop[tr]['votemper_m'] - gsw.CT_freezing(prop[tr]['vosaline_m'],depth, air_frac) - prop[tr]['votemper_e'] +gsw.CT_freezing(prop[tr]['vosaline_e'],depth, air_frac) + prop[tr]['votemper_std_e']
                                        dTf_me_i[t,COORDy[data[tr]['eddy_j']].astype(int),COORDx[data[tr]['eddy_i']].astype(int)] = prop[tr]['votemper_m'] - gsw.CT_freezing(prop[tr]['vosaline_m'],depth, air_frac) - prop[tr]['votemper_e'] +gsw.CT_freezing(prop[tr]['vosaline_e'],depth, air_frac) + prop[tr]['votemper_std_e'] + gsw.CT_freezing(prop[tr]['vosaline_e']-prop[tr]['vosaline_std_e'],depth, air_frac)-gsw.CT_freezing(prop[tr]['vosaline_e']+prop[tr]['vosaline_std_e'],depth, air_frac)
                    else:
                        print('No detection on day ', timeframe)

                dstemp = de.copy(deep = True)
                dstemp = dstemp.squeeze()
                dstemp['day'] = (('day'), np.arange(nmb_days))
                dstemp = dstemp.set_coords('day')
                dstemp.drop('time_counter')

                dstemp['Nb'] = (('day', 'y', 'x'), Nb)
                dstemp['T_m_i'] = (('day', 'y', 'x'), T_m_i)
                dstemp['S_m_i'] = (('day', 'y', 'x'), S_m_i)
                dstemp['DT_me_i'] = (('day', 'y', 'x'), DT_me_i)
                dstemp['DS_me_i'] = (('day', 'y', 'x'), DS_me_i)
                dstemp['dT_me_i'] = (('day', 'y', 'x'), dT_me_i)
                dstemp['dS_me_i'] = (('day', 'y', 'x'), dS_me_i)
                dstemp['dTf_me_i'] = (('day', 'y', 'x'), dTf_me_i)

                dstemp['dtb'] = (('day', 'y', 'x'), dtb)


                dsAll['Nb_sum'] = dstemp.Nb.where(dstemp.dtb>=filter_time).sum('day', skipna = True).astype(np.uint16).expand_dims('time_counter')

                dsAll['T_m_i'] = dstemp.T_m_i.where(dstemp.dtb>=filter_time).mean('day', skipna = True).expand_dims('time_counter')

                dsAll['S_m_i'] = dstemp.S_m_i.where(dstemp.dtb>=filter_time).mean('day', skipna = True).expand_dims( 'time_counter')
                
                dsAll['DT_me_i'] = dstemp.DT_me_i.where(dstemp.dtb>=filter_time).mean('day', skipna = True).expand_dims( 'time_counter')

                dsAll['DS_me_i'] = dstemp.DS_me_i.where(dstemp.dtb>=filter_time).mean('day', skipna = True).expand_dims( 'time_counter')

                dsAll['dT_me_i'] = dstemp.dT_me_i.where(dstemp.dT_me_i*dstemp.DT_me_i>0).where(dstemp.dtb>=filter_time).mean('day', skipna = True).expand_dims( 'time_counter')

                dsAll['dS_me_i'] = dstemp.dS_me_i.where(dstemp.dS_me_i*dstemp.DS_me_i>0).where(dstemp.dtb>=filter_time).mean('day', skipna = True).expand_dims( 'time_counter')
                
                dsAll['dTf_me_i'] = dstemp.dTf_me_i.where(dstemp.dT_me_i*dstemp.DT_me_i>0).where(dstemp.dtb>=filter_time).mean('day', skipna = True).expand_dims( 'time_counter')

                #save ds0 y-m
                dsAll.to_netcdf(Config['save_folder_path']+Config['prefixe_monthyP_means']+str(y)+'_'+'{0:02}'.format(m)+'_depthi'+str(depth_ind)+'_'+Config['suffixe_all']+'.nc')
            else:
                print('I already have dsAll for y: ', y, ' m : ', m)
    return None

