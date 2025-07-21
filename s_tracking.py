import sys
import f_tracking


depth_index = int(eval(sys.argv[1]))
from Config_file_REF12_3 import Config, Model  

yi = 1995 
yf = 2020

print('I start tracking')
f_tracking.tracking_REF12(depth_index,yi, yf, Config)
f_tracking.get_eddies_1d(depth_index, yi, yf, Config)