import sys
import f_interp


depth_index = int(eval(sys.argv[1]))
yi = int(eval(sys.argv[2]))
str_conf = sys.argv[3]

if str_conf =='G':
    from Config_file_G12 import Config, Model 
    convert = False 
elif str_conf == 'C':
    from Config_file_REF12_3 import Config, Model  
    convert = False
 
css = Config['css'] # type: ignore
f_interp.interpolate_new_rolling(depth_index, yi, yi, Config, Model, lsig = 50, cs_vert = css[0], cs_x = css[1], cs_y = css[2], should_convert_cftime=convert) # type: ignore