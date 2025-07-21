import cftime
import xarray as xr
from glob import glob
import numpy as np
import os
import pickle
from xorca.lib import load_xorca_dataset # type: ignore
import xgcm # type: ignore 
import f_generic
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.simplefilter('ignore', UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)


def Prop_REF12_TSonly_m(i_depth, Config, y, m):
    fact = Config['properties_parameters']['fact']   
    ds_zgr = xr.open_dataset(Config['domain_path'])
    L_depth_value = ds_zgr.nav_lev.values.astype(int)

    if not os.path.exists(Config['save_folder_path'] + Config['prefixe_properties_m'] +str(y)+'_'+'{0:02}'.format(m)+ '_depthi' + str(i_depth) + Config['suffixe_all']+'.pickle'):

        Liste_mesh = [Config['domain_path']]
        dsT = load_xorca_dataset(data_files=[Config['datapath']+'CREG12.L75-REF12_y'+str(y)+'m'+'{0:02}'.format(m)+'.1d_gridT.nc'],model_config='NEST', aux_files=Liste_mesh)
        dsT['at'] = dsT['e1t']*dsT['e2t']
        dsT['au'] = dsT['e1u']*dsT['e2u']
        dsT['av'] = dsT['e1v']*dsT['e2v']
        dsT['af'] = dsT['e1f']*dsT['e2f']
        dsT = dsT.set_coords(('at', 'au', 'av', 'af'))
        dsT = dsT.isel(y_r = slice(0, -1)).isel(y_c = slice(0, -1)).isel(x_c=slice(None, -2)).isel(x_r=slice(None, -2)).isel(z_c = i_depth)


    
        metrics = {('X',): ['e1u', 'e1v'], ('Y',): ['e2v', 'e2u'], ('X', 'Y'): ['at', 'au', 'av', 'af']}
        coords={"X": {"center": "x_c", "right": "x_r"}, "Y":{'center':'y_c', 'right':'y_r'}}
        ds_grid = f_generic.save_grid_REF12(Config)

        def get_ind_env(eddy_i_ind, eddy_j_ind, fact, j_max = dsT.y_r.max().values, i_max = dsT.x_r.max().values):
            eddy_i_size = int((eddy_i_ind.max() - eddy_i_ind.min()).values/2)
            eddy_j_size = int((eddy_j_ind.max() - eddy_j_ind.min()).values/2)
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

        gridT = xgcm.Grid(dsT, periodic=["X"], coords = coords, metrics=metrics, boundary = 'fill', fill_value=np.nan)

        # interpolate fields on f-points
        print('Interpolate fields on f-points')
        votemper_f = gridT.interp(dsT.votemper, axis = ('X', 'Y')).compute()
        vosaline_f = gridT.interp(dsT.vosaline, axis = ('X', 'Y')).compute()

        print('Get values per eddy per day')
        Properties_year = {}
        file= Config['save_folder_path'] +Config['prefixe_detection_m'] +str(y)+'_'+'{0:02}'.format(m)+'_depthi'+str(i_depth)+Config['suffixe_all']+'_eddies.pickle'
 
        infile = open(file,'rb')
        detected = pickle.load(infile)
        infile.close()


        for d in range(len(detected)):
            print(d, '/', dsT.t.shape[0])
            properties_per_day = {}
            frame = detected[d]
            if len(frame)>=1:
                temp_d = votemper_f.sel(t =frame[0]['time'])
                sal_d = vosaline_f.sel(t=frame[0]['time'])


                for e in range(len(frame)):
                    eddy_j_ind = ds_grid.y[frame[e]['eddy_j']]
                    eddy_i_ind = ds_grid.x[frame[e]['eddy_i']]
                    # function
                    temp = temp_d.sel(x_r = eddy_i_ind).sel(y_r = eddy_j_ind).mean()
                    sal = sal_d.sel(x_r = eddy_i_ind).sel(y_r = eddy_j_ind).mean()                    
                    #anomaly
                    i_env, j_env = get_ind_env(eddy_i_ind, eddy_j_ind, fact) 

                    #temp
                    T_centre = temp_d.sel(x_r =sorted(eddy_i_ind.values)[len(eddy_i_ind)//2]).sel(y_r = sorted(eddy_j_ind.values)[len(eddy_j_ind)//2]).mean()
                    T_env = temp_d.sel(x_r = i_env).sel(y_r = j_env).mean()
                    T_std_env = temp_d.sel(x_r = i_env).sel(y_r = j_env).std()
                    #sal
                    S_centre = sal_d.sel(x_r =sorted(eddy_i_ind.values)[len(eddy_i_ind)//2]).sel(y_r = sorted(eddy_j_ind.values)[len(eddy_j_ind)//2]).mean()
                    S_env = sal_d.sel(x_r = i_env).sel(y_r = j_env).mean()
                    S_std_env = sal_d.sel(x_r = i_env).sel(y_r = j_env).std()
                    # Write all 
                    properties_per_day[e] = {'time':frame[0]['time'],
                                            'votemper_m':temp.values, 'vosaline_m':sal.values,
                                            'votemper_c':T_centre.values, 'votemper_e':T_env.values,
                                            'vosaline_c':S_centre.values, 'vosaline_e':S_env.values,
                                            'vosaline_std_e':S_std_env, 'votemper_std_e':T_std_env,
}


                Properties_year[d] = properties_per_day
        # save
        with open(Config['save_folder_path'] + Config['prefixe_properties_m'] +str(y)+'_'+'{0:02}'.format(m)+ '_depthi' + str(i_depth) + Config['suffixe_all']+'.pickle', 'wb') as f:
            pickle.dump(Properties_year, f, pickle.HIGHEST_PROTOCOL)
        f.close()
    else:
        print('I already have properties for ', y, '/', m)
    return None



def Prop_REF12_Iceonly_m(i_depth, Config, y, m):
    fact = Config['properties_parameters']['fact']   
    ds_zgr = xr.open_dataset(Config['domain_path'])
    L_depth_value = ds_zgr.nav_lev.values.astype(int)

    if not os.path.exists(Config['save_folder_path'] + Config['prefixe_properties_m'] +str(y)+'_'+'{0:02}'.format(m)+ '_depthi' + str(i_depth) + Config['suffixe_all']+'.pickle'):

        Liste_mesh = [Config['domain_path']]
        dsI = load_xorca_dataset(data_files=[Config['datapath']+'CREG12.L75-REF12_y'+str(y)+'m'+'{0:02}'.format(m)+'.1d_icemod.nc'],
                                update_orca_variables = {'sivolu' : {"dims":["t", "y_c", "x_c",]}, 
                                                        'siconc' : {"dims":["t", "y_c", "x_c",]},
                                                        },
                                model_config='NEST', aux_files=Liste_mesh)
        dsI['at'] = dsI['e1t']*dsI['e2t']
        dsI['au'] = dsI['e1u']*dsI['e2u']
        dsI['av'] = dsI['e1v']*dsI['e2v']
        dsI['af'] = dsI['e1f']*dsI['e2f']
        dsI = dsI.set_coords(('at', 'au', 'av', 'af'))
        dsI = dsI.isel(y_r = slice(0, -1)).isel(y_c = slice(0, -1)).isel(x_c=slice(None, -2)).isel(x_r=slice(None, -2))

        metrics = {('X',): ['e1u', 'e1v'], ('Y',): ['e2v', 'e2u'], ('X', 'Y'): ['at', 'au', 'av', 'af']}
        coords={"X": {"center": "x_c", "right": "x_r"}, "Y":{'center':'y_c', 'right':'y_r'}}
        ds_grid = f_generic.save_grid_REF12(Config)

        def get_ind_env(eddy_i_ind, eddy_j_ind, fact, j_max = dsI.y_r.max().values, i_max = dsI.x_r.max().values):
            eddy_i_size = int((eddy_i_ind.max() - eddy_i_ind.min()).values/2)
            eddy_j_size = int((eddy_j_ind.max() - eddy_j_ind.min()).values/2)
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

        gridI = xgcm.Grid(dsI, periodic=["X"], coords = coords, metrics=metrics, boundary = 'fill', fill_value=np.nan)

        # interpolate fields on f-points
        print('Interpolate fields on f-points')
        sivolu_f = gridI.interp(dsI.sivolu, axis = ('X', 'Y')).compute()
        siconc_f = gridI.interp(dsI.siconc, axis = ('X', 'Y')).compute()

        print('Get values per eddy per day')
        Properties_year = {}
        file= Config['save_folder_path'] +Config['prefixe_detection_m'] +str(y)+'_'+'{0:02}'.format(m)+'_depthi'+str(i_depth)+Config['suffixe_all']+'_eddies.pickle'
 
        infile = open(file,'rb')
        detected = pickle.load(infile)
        infile.close()


        for d in range(len(detected)):
            print(d, '/', dsI.t.shape[0])
            properties_per_day = {}
            frame = detected[d]
            if len(frame)>=1:
                fra_d = siconc_f.sel(t=frame[0]['time'])
                thic_d = sivolu_f.sel(t=frame[0]['time'])


                for e in range(len(frame)):
                    eddy_j_ind = ds_grid.y[frame[e]['eddy_j']]
                    eddy_i_ind = ds_grid.x[frame[e]['eddy_i']]
                    # function
                    fra = fra_d.sel(x_r = eddy_i_ind).sel(y_r = eddy_j_ind).mean()
                    thic = thic_d.sel(x_r = eddy_i_ind).sel(y_r = eddy_j_ind).mean()
               
                    #anomaly
                    i_env, j_env = get_ind_env(eddy_i_ind, eddy_j_ind, fact) 

                    #Ice conc
                    If_env = fra_d.sel(x_r = i_env).sel(y_r = j_env).mean()
                    Ih_env = thic_d.sel(x_r = i_env).sel(y_r = j_env).mean()
                    If_std = fra_d.sel(x_r = i_env).sel(y_r = j_env).std()
                    Ih_std = thic_d.sel(x_r = i_env).sel(y_r = j_env).std()

                    # Write all 
                    properties_per_day[e] = {'time':frame[0]['time'],
                                            'ifra_m':fra.values, 'ithic_m':thic.values,
                                            'ifra_e':If_env.values, 'ithic_e':Ih_env.values,
                                            'ifra_std':If_std.values, 'ithic_std':Ih_std.values,
}


                Properties_year[d] = properties_per_day
        # save
        with open(Config['save_folder_path'] + Config['prefixe_properties_m'] +str(y)+'_'+'{0:02}'.format(m)+ '_depthi' + str(i_depth) + Config['suffixe_all']+'.pickle', 'wb') as f:
            pickle.dump(Properties_year, f, pickle.HIGHEST_PROTOCOL)
        f.close()
    else:
        print('I already have properties for ', y, '/', m)
    return None

