import numpy as np
# version restricted to the CB 

Model = {'model':'ORCA',
        'grid':'latlon',
        'calendar':'NoLeap',
        'lon1': 180, # minimum longitude of detection region /!\ attention Harcode dans Interp
        'lon2': -180, # maximum longitude /!\ attention Harcode dans Interp
        'lat1': 67, # minimum latitude /!\ attention Harcode dans Interp
        'lat2': 90, # maximum latitude /!\ attention Harcode dans Interp
        'res' :  1./12., # resolution of the fields in degrees
        'i1': 200, # min lon in grid coordinates
        'i2' :1150,# max lon in grid coordinates
        'j1':1250,# min lat in grid coordinates
        'j2':1700 # max lat in grid coordinates
}

Config  = {'suffixe_all':'REF12_3',
           'L_depth_index':np.arange(0,75), 
            'domain_path':'/home/nplanat/projects/ctb-cdufour/nplanat/CREG_mask_modifs/CREG12.L75-REF09_domain_cfg_20230801_Z.nc',
            'prefixe_interpolation':'Interp/Interp_CREG_1d_std50_',
            'prefixe_interpolation_c':'OW',
            'prefixe_detection_m':'Detected/Detected_CREG_1d_OWthr_03_m_',
            'prefixe_tracking':'Tracked/Tracked_CREG_1d_1995_2020_',
            'prefixe_tracking_det_indexes_m':'Track_Det/Track_Det_Index_CREG_1d_1995_2020_m_',
            'tracking_formated_prefixe':'Formatted/1995_2020_Ice',
            'prefixe_monthy_means':'Monthly/M_',
            'prefixe_monthyP_means':'Monthly/PIce_',
            'prefixe_monthly_climato':'Monthly/C_',
            'prefixe_properties_m' : 'Prop/Prop_ice_',
            'prefixe_properties_tracking':'Prop_Tracked/Tracked_CREG_1d_1995_2020_Ice_',
            'save_folder_path':'/home/nplanat/scratch/REF12_3/',
            'datapath':'/home/nplanat/projects/ctb-cdufour/shared_files/models/CREG12-REF12/1d/', 
            'datapath_c':'/home/nplanat/projects/ctb-cdufour/shared_files/models/CREG12-REF12/climato/', 
            'meshpath':'/home/nplanat/projects/ctb-cdufour/nplanat/CREG_mask_modifs/CREG12.L75-REF09_domain_cfg_20230801_Z.nc',
            'maskpath':'/home/nplanat/projects/ctb-cdufour/nplanat/CREG_mask_modifs/CREG_correction_latlon_mask.nc',
            'one_T_file_path' : '/home/nplanat/projects/ctb-cdufour/shared_files/models/CREG12-REF12/1d/CREG12.L75-REF12_y2003m03.1d_gridT.nc',
            'interpolation_parameters':{'model': Model['model'],
                                        'grid': Model['grid'],
                                        'calendar': Model['calendar'], # calendar, must be either 360_day or standard or NoLeap
                                        'lon1': Model['lon1'], # minimum longitude of detection region
                                        'lon2': Model['lon2'], # maximum longitude
                                        'lat1': Model['lat1'], # minimum latitude
                                        'lat2': Model['lat2'], # maximum latitude
                                        'res': Model['res'], # resolution of the fields in degrees
                                        'vars_to_interpolate': ['OW', 'vort'], # variables to be interpolated 
                                        'mask_to_interpolate': ['fmask', 'tmask', 'bathymetry'], # masks to interpolate
                                        'regrid_method': 'bilinear', # method used for regridding (default is 'bilinear')
                                        'ext_method': None}, # masks to interpolate\
           'detection_parameters':{'model': Model['model'],
                                    'grid': Model['grid'],
                                    'calendar': Model['calendar'], # calendar, must be either 360_day or standard
                                    'lon1': Model['lon1'], # minimum longitude of detection region
                                    'lon2': Model['lon2'], # maximum longitude
                                    'lat1': Model['lat1'], # minimum latitude
                                    'lat2': Model['lat2'], # maximum latitude
                                    'res': Model['res'], # resolution of the fields in degrees
                                    'OW_thr_name': 'OW_std', # Okubo-Weiss threshold for eddy detection
                                    'OW_thr_factor': -0.3, # Okubo-Weiss parameter threshold
                                    'Npix_min': 20, # minimum number of pixels (grid cells) to be considered as eddy
                                    'Npix_max': 500,
                                    'no_long': False, # If True, elongated shapes will not be considered
                                    'no_two': False, # If True, eddies with two minima in the OW
                                                        # parameter and a OW > OW_thr in between  will not
                                                        # be considered} # maximum number of pixels (grid cells)
                                    'suffixe_OW':'.1995_2020',
                                    },
            'tracking_parameters':{'model': Model['model'],
                                    'grid': Model['grid'],
                                    'dt': 1,
                                    'calendar' : Model['calendar'], # calendar, must be either 360_day or standard or NoLeap
                                    'dates_of_detection':None, 
                                    'lon1': Model['lon1'], # minimum longitude of detection region
                                    'lon2': Model['lon2'], # maximum longitude
                                    'lat1': Model['lat1'], # minimum latitude
                                    'lat2': Model['lat2'], # maximum latitude
                                    'search_dist': 22000, # maximum distance of search ellipsis from eddy center in towards the east 
                                             # (if set to 0, it will be calculated as (150. / (7. / dt)))
                                   'search_circle': True,
                                   'eddy_scale_min': 0.5, # minimum factor by which eddy amplitude and area are allowed to change in one timestep
                                   'eddy_scale_max': 1.5, # maximum factor by which eddy amplitude and area are allowed to change in one timestep
                                   'dict': 0, # dictionary containing detected eddies to be used when not stored in files (set to 0 otherwise)
                                   'ross_path': '',
                                   'detection_is_saved_monthly_files':True #, True os we save monthly files for detection - which I do
                                   },
            'filtering_parameters':{0:{'name':'Amerasian_dt2_bathy0',
                                    'lon1':-108, # keep eddies <= lon1
                                    'lon2':180, # kepp eddies >= lon2
                                    'lat1':69, # keep eddies >= lat1
                                    'lat2':85, # keep eddies <= lat2
                                    'bathy_min':0,# keep eddies where bathy >= bathy_min
                                    'yi':1995, #decade init
                                    'yf':2020, #decade fin
                                    'duration':2 #number of days minimum
                                    }},
            'properties_parameters':{'fact':3, 
                                    'Properties_extract' : ['ifra_m', 'ithic_m',\
                                                            'ifra_e', 'ithic_e',\
                                                            'ifra_std', 'ithic_std']},
            'forcings_parameters':{'w_PV_background':61,
                                   'XMARGINS':21},
}