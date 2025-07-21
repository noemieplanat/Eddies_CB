
import xarray as xr
import os
import numpy as np
import gc

def format_monthly_means_glouton_simplified(ind, Config, yii, yff):
    ds_zgr = xr.open_dataset(Config['domain_path'])
    L_depth_value = ds_zgr.nav_lev.values.astype(int)
    Nav_Lon = ds_zgr.nav_lon.where(ds_zgr.nav_lon>=0, ds_zgr.nav_lon+360)
    ds_zgr['nav_lon360'] = Nav_Lon
    ds_zgr = ds_zgr.set_coords({'nav_lon', 'nav_lat', 'nav_lon360'})
    
    if not os.path.exists(Config['save_folder_path']+Config['prefixe_monthly_climato']+str(yii)+'_'+str(yff)+'depthi_'+str(ind)+Config['suffixe_all']+'.nc'):
        ds0 = load_monthly_mean_fields(ind, Config,yii, yff)
        print(ds0.time_counter.shape[0],  (yff+1-yii)*12)
        if ds0.time_counter.shape[0] == (yff+1-yii)*12:
            ds = ds0.where(ds0.time_counter.dt.year>= yii, drop = True).where(ds0.time_counter.dt.year<=yff, drop = True)
            dss = ds.groupby(ds.time_counter.dt.month).sum(dim = 'time_counter')
            ds1 = ds.copy(deep = True)
            Dec1 = ds1.groupby(ds1.time_counter.dt.month).mean(dim = 'time_counter')
            Dec1['Nb_sum'] = dss.Nb_sum
            Dec1['NCb_sum'] = dss.NCb_sum
            ds2 = ds.copy(deep = True)
            ds3 = ds.copy(deep = True)
            ds3 = ds3.where(ds3.Nb_sum>0, np.nan)
            Dec = ds3.groupby(ds3.time_counter.dt.month).mean('time_counter')
            Dec1['Rb_mean'] = Dec['Rb_mean']
            Dec1['Ab_mean'] = Dec['Ab_mean']
            Dec1['dtb_mean'] = Dec['dtb_mean']
            Dec1['dtb_mean'] = Dec['dtb_mean']
            Dec1['D'] = Dec['D']

            Dec1['nav_lon'] = (('y', 'x'), ds0.nav_lon.isel(time_counter = 0).values)
            Dec1['nav_lat'] = (('y', 'x'), ds0.nav_lat.isel(time_counter = 0).values)
            #Dec1['nav_lev'] = ds_zgr['nav_lev'] 
            Dec1.to_netcdf(Config['save_folder_path']+Config['prefixe_monthly_climato']+str(yii)+'_'+str(yff)+'depthi_'+str(ind)+Config['suffixe_all']+'.nc')
            del Dec1
        else:
            print('I dont have all monthly files for ind ', ind)
    return None  

def format_monthly_prop_glouton_simplified(ind, Config, yii, yff):
    ds_zgr = xr.open_dataset(Config['domain_path'])
    L_depth_value = ds_zgr.nav_lev.values.astype(int)
    Nav_Lon = ds_zgr.nav_lon.where(ds_zgr.nav_lon>=0, ds_zgr.nav_lon+360)
    ds_zgr['nav_lon360'] = Nav_Lon
    ds_zgr = ds_zgr.set_coords({'nav_lon', 'nav_lat', 'nav_lon360'})
    
    if not os.path.exists(Config['save_folder_path']+Config['prefixe_monthlyP_climato']+str(yii)+'_'+str(yff)+'depthi_'+str(ind)+Config['suffixe_all']+'.nc'):
        ds0 = load_monthly_mean_prop(ind, Config,yii, yff)
        print(ds0.time_counter.shape[0],  (yff+1-yii)*12)
        if ds0.time_counter.shape[0] == (yff+1-yii)*12:
            ds = ds0.where(ds0.time_counter.dt.year>= yii, drop = True).where(ds0.time_counter.dt.year<=yff, drop = True)
            dss = ds.groupby(ds.time_counter.dt.month).sum(dim = 'time_counter')
            ds1 = ds.copy(deep = True)
            Dec1 = ds1.groupby(ds1.time_counter.dt.month).mean(dim = 'time_counter')
            Dec1['Nb_sum'] = dss.Nb_sum
            ds3 = ds.copy(deep = True)
            ds3 = ds3.where(ds3.Nb_sum>0, np.nan)
            Dec = ds3.groupby(ds3.time_counter.dt.month).mean('time_counter')
            Dec1['T_m_i'] = Dec['T_m_i'] # erreur : Rb_mean
            Dec1['S_m_i'] = Dec['S_m_i'] # erreur : Ab_mean
            Dec1['DT_me_i'] = Dec['DT_me_i']
            Dec1['DS_me_i'] = Dec['DS_me_i']
            Dec1['dT_me_i'] = Dec['dT_me_i']
            Dec1['dS_me_i'] = Dec['dS_me_i']
            Dec1['dTf_me_i'] = Dec['dTf_me_i']

            Dec1['nav_lon'] = (('y', 'x'), ds0.nav_lon.isel(time_counter = 0).values)
            Dec1['nav_lat'] = (('y', 'x'), ds0.nav_lat.isel(time_counter = 0).values)
            #Dec1['nav_lev'] = ds_zgr['nav_lev'] 
            Dec1.to_netcdf(Config['save_folder_path']+Config['prefixe_monthlyP_climato']+str(yii)+'_'+str(yff)+'depthi_'+str(ind)+Config['suffixe_all']+'.nc')
            del Dec1
        else:
            print('I dont have all monthly files for ind ', ind)
    return None 

def load_monthly_mean_fields(ind, Config,yi, yf):
    Liste_monthly_files = []
    for y in range(yi, yf+1):
        for month in range(1, 13):
            Liste_monthly_files.append(Config['save_folder_path']+Config['prefixe_monthy_means']+str(y)+'_'+'{0:02}'.format(month)+'_depthi'+str(ind)+'_'+Config['suffixe_all']+'.nc')                
    dsAll = xr.open_dataset(Liste_monthly_files[0])
    dsAll = dsAll.set_coords({'nav_lat', 'nav_lon'})
    for i in range(1, len(Liste_monthly_files)):
        print(i, '/', len(Liste_monthly_files))
        ds = xr.open_dataset(Liste_monthly_files[i])
        ds = ds.set_coords({'nav_lat', 'nav_lon'})
        dsAll = xr.concat((dsAll, ds), dim = 'time_counter')
    return dsAll

def load_monthly_mean_prop(ind, Config,yi, yf):
    Liste_monthly_files = []
    for y in range(yi, yf+1):
        for month in range(1, 13):
            Liste_monthly_files.append(Config['save_folder_path']+Config['prefixe_monthyP_means']+str(y)+'_'+'{0:02}'.format(month)+'_depthi'+str(ind)+'_'+Config['suffixe_all']+'.nc')                
    dsAll = xr.open_dataset(Liste_monthly_files[0])
    dsAll = dsAll.set_coords({'nav_lat', 'nav_lon'})
    for i in range(1, len(Liste_monthly_files)):
        print(i, '/', len(Liste_monthly_files))
        ds = xr.open_dataset(Liste_monthly_files[i])
        ds = ds.set_coords({'nav_lat', 'nav_lon'})
        dsAll = xr.concat((dsAll, ds), dim = 'time_counter')
    return dsAll