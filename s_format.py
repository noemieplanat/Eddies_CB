import sys
import f_format


depth_index = int(eval(sys.argv[1]))
filter_format = int(eval(sys.argv[2]))


from Config_file_REF12_3_short_track import Config, Model  
    
f_format.format_save_tracked_filtered(Config, depth_index, filter_format)
f_format.format_save_tracked_prop_filtered(Config, depth_index, filter_format)
