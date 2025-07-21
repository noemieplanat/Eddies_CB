import sys
import f_detect

yi = int(eval(sys.argv[2]))
depth_index = int(eval(sys.argv[1]))

from Config_file_REF12_3 import Config, Model  
f_detect.detection_save_m(depth_index, yi, yi, Config, 12, Model) #yf is included ! 