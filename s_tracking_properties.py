import sys
import f_tracking


depth_index = int(eval(sys.argv[1]))

from Config_file_REF12_3 import Config, Model  
yi = 1995; yf = 2020

f_tracking.track_properties_REF12(depth_index, Config, yi, yf)
f_tracking.track_properties_REF12_1d(depth_index,Config, yi, yf)