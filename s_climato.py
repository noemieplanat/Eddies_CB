import sys
import f_climato

depth_index = int(eval(sys.argv[1]))

from Config_file_REF12_3 import Config, Model  

yi = 1995; yf = 2020
f_climato.format_monthly_means_glouton_simplified(depth_index, Config, yi, yf)
f_climato.format_monthly_prop_glouton_simplified(depth_index, Config, yi, yf)
