from glob import glob
import xarray as xr
import xgcm
import numpy as np
import pickle
import os
from xorca.lib import load_xorca_dataset
import f_generic
import cftime as cf

def get_gridE(y,m,Config):
    if not os.path.exists(Config['datapath']+'/CREG12.L75-REF12_y'+str(y)+'m{0:02}'.format(m)+'.1d_gridE.nc'):
        #load a file with xorca to create a grid object and interpolate
        one_T_file_path=Config['datapath']+'/CREG12.L75-REF12_y'+str(y)+'m{0:02}'.format(m)+'.1d_gridT.nc'
        domain_path = Config['domain_path']
        Liste_mesh = [domain_path]
        dsT = load_xorca_dataset(data_files=[one_T_file_path],model_config='NEST', aux_files=Liste_mesh)
        dsT['at'] = dsT['e1t']*dsT['e2t']
        dsT['au'] = dsT['e1u']*dsT['e2u']
        dsT['av'] = dsT['e1v']*dsT['e2v']
        dsT['af'] = dsT['e1f']*dsT['e2f']
        dsT = dsT.set_coords(('at', 'au', 'av', 'af'))
        dsT = dsT.isel(y_r = slice(0, -1)).isel(y_c = slice(0, -1)).isel(x_c=slice(None, -2)).isel(x_r=slice(None, -2))

        metrics = {('X',): ['e1u', 'e1v'], ('Y',): ['e2v', 'e2u'], ('X', 'Y'): ['at', 'au', 'av', 'af']}
        coords={"X": {"center": "x_c", "right": "x_r"}, "Y":{'center':'y_c', 'right':'y_r'}}


        ds_grid = f_generic.save_grid_REF12(Config)

        COORDx = np.zeros((dsT.x_r.shape[0], ))
        for i in range(ds_grid.x.shape[0]):
            ind = np.where(dsT.x_r.values == ds_grid.x.values[i])[0]
            COORDx[i] = ind

        COORDy = np.zeros((dsT.y_r.shape[0], ))
        for i in range(ds_grid.y.shape[0]):
            ind = np.where(dsT.y_r.values == ds_grid.y.values[i])[0]
            COORDy[i] = ind
        gridT = xgcm.Grid(dsT, periodic=["X"], coords = coords, metrics=metrics, boundary = 'fill', fill_value=np.nan)


        dsM = xr.open_dataset(Config['Gridf'])
        if m in [1,3,5,7,8,10,12]:
            dm = 31
        elif m in [4,6,9,11]:
            dm = 30
        else:
            dm = 28
        dsM['time_counter'] = (('time_counter'), [cf.datetime(y,m,d,12, calendar = 'noleap') for d in range(1, dm+1)])
        Mask = np.zeros((dsT.t.shape[0],dsT.z_c.shape[0], dsT.y_c.shape[0]+1, dsT.x_c.shape[0]+2), dtype = np.int8) 
        R = np.zeros((dsT.t.shape[0],dsT.z_c.shape[0], dsT.y_c.shape[0]+1, dsT.x_c.shape[0]+2), dtype = 'f')  
        A = np.zeros((dsT.t.shape[0],dsT.z_c.shape[0], dsT.y_c.shape[0]+1, dsT.x_c.shape[0]+2), dtype = 'f')  
        print('get eddies')
        for depth_index in range(54):
            if os.path.exists(Config['save_folder_path'] +Config['prefixe_detection_m']+str(y)+'_{0:02}'.format(m)+'_depthi'+str(depth_index)+Config['suffixe_all']+'_eddies.pickle'):
                #try:
                infile = open(Config['save_folder_path'] +Config['prefixe_detection_m']+str(y)+'_{0:02}'.format(m)+'_depthi'+str(depth_index)+Config['suffixe_all']+'_eddies.pickle','rb')
                data = pickle.load(infile)
                infile.close()
                for it in range(dsM.time_counter.shape[0]):
                    for e in range(len(data[it])):
                        if data[it][e]['type'] == 'cyclonic':
                            Mask[it,depth_index,1+COORDy[data[it][e]['eddy_j']].astype(int),1+COORDx[data[it][e]['eddy_i']].astype(int)] = -1
                        else:
                            Mask[it,depth_index,1+COORDy[data[it][e]['eddy_j']].astype(int),1+COORDx[data[it][e]['eddy_i']].astype(int)] = 1
                        R[it,depth_index,1+COORDy[data[it][e]['eddy_j']].astype(int),1+COORDx[data[it][e]['eddy_i']].astype(int)] = data[it][e]['scale']
                        A[it,depth_index,1+COORDy[data[it][e]['eddy_j']].astype(int),1+COORDx[data[it][e]['eddy_i']].astype(int)] = data[it][e]['amp']
                #except:
                #    print('Problem with time frame ', str(y)+'-{0:02}'.format(m), 'at ind depth ', depth_index)
            else:
                print('No eddy detected at ind depth ', depth_index, 'on month ', m)

        dsM['Mask'] = (('time_counter', 'deptht', 'y', 'x'),Mask)
        del Mask
        print('R')
        dsM['R'] = (('time_counter', 'deptht', 'y', 'x'),R)
        del R
        print('A')
        dsM['A'] = (('time_counter', 'deptht', 'y', 'x'),A)
        del A
        print('saving')
        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in dsM.data_vars}
        dsM.to_netcdf(Config['datapath']+'/CREG12.L75-REF12_y'+str(y)+'m{0:02}'.format(m)+'.1d_gridE.nc', encoding = encoding)
    else:
        print('I already have gridE for ', str(y)+'m{0:02}'.format(m))
        