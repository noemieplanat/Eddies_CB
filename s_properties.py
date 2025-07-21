import sys
import f_properties


y = int(eval(sys.argv[1]))
depth_index = int(eval(sys.argv[2]))
m = int(eval(sys.argv[3]))

from Config_file_REF12_3_prop_ice import Config, Model  
    
print('Doing ', m, '/', y)
f_properties.Prop_REF12_Iceonly_m(depth_index, Config, y, m)