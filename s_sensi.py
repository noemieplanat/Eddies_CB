import sys
import f_detect
import f_interp
import f_tracking
import f_format

sensi = int(eval(sys.argv[1]))


idepths0 = [14]#, 20,27,32,39]
y0 = 2000

from Config_file_REF12_3 import Config, Model

if sensi ==1:
    Config['detection_parameters']['OW_thr_factor'] = -0.1
    lsig = 50
if sensi ==2:
    Config['detection_parameters']['OW_thr_factor'] = -0.2
    lsig = 50
if sensi ==3:
    Config['detection_parameters']['OW_thr_factor'] = -0.3
    lsig = 50
if sensi ==4:
    Config['detection_parameters']['OW_thr_factor'] = -0.4
    lsig = 50
if sensi ==5:
    Config['detection_parameters']['OW_thr_factor'] = -0.5
    lsig = 50
if sensi ==6:
    Config['detection_parameters']['OW_thr_factor'] = -0.3
    lsig = 24
    Config['detection_parameters']['suffixe_OW'] = '.2000'
if sensi ==7:
    Config['detection_parameters']['OW_thr_factor'] = -0.3
    lsig = 74
    Config['detection_parameters']['suffixe_OW'] = '.2000'
if sensi ==8:
    Config['detection_parameters']['OW_thr_factor'] = -0.3
    lsig = 100
    Config['detection_parameters']['suffixe_OW'] = '.2000'
if sensi ==9:
    Config['detection_parameters']['OW_thr_factor'] = -0.3
    lsig = 12
    Config['detection_parameters']['suffixe_OW'] = '.2000'

Config['prefixe_detection_m']= 'Detected/Detected_CREG_1d_OWthr_SENSI_'+str(Config['detection_parameters']['OW_thr_factor'])[1:]+'_m_'+str(lsig)+'_'
Config['prefixe_interpolation'] = 'Interp/Interp_CREG_1d_std'+str(lsig)+'_'
Config['prefixe_tracking'] = 'Tracked/Tracked_CREG_1d_SENSI_'+str(sensi)
Config['prefixe_tracking_det_indexes_m']= 'Track_Det/Track_Det_Index_CREG_1d_SENSI_'+str(sensi)
Config['tracking_formated_prefixe'] = 'Formatted/SENSI_'+str(sensi)




for depth_index in idepths0:
    print(depth_index, 'interp')
    if sensi>5:
        f_interp.interpolate_new_rolling(depth_index, y0, y0, Config, Model, lsig = lsig)
        print(depth_index, 'detec')
        f_detect.detection_save_m(depth_index, y0, y0, Config, 12, Model, old_interp = False) #yf is included ! 
    else:
        print(depth_index, 'track')
        f_detect.detection_save_m(depth_index, y0, y0, Config, 12, Model, old_interp = True) #yf is included ! 
    f_tracking.tracking_REF12(depth_index,y0, y0, Config)
    print(depth_index, '1d')
    f_tracking.get_eddies_1d(depth_index, y0, y0, Config)
    print(depth_index, 'format filtered')
    f_format.format_save_tracked_filtered(Config, depth_index, 0)


