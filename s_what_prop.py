from Config_file_REF12_3 import Config
import os

for idepth in range(0, 54):
    have_all_depth = True
    for y in range(1995, 2021):
        for m in range(1, 13):
            if not os.path.exists(Config['save_folder_path']+Config['prefixe_properties_m']+str(y)+'_{0:02}'.format(m)+ '_depthi' + str(idepth) + Config['suffixe_all']+'.pickle')\
                  or os.path.getsize(Config['save_folder_path']+Config['prefixe_properties_m']+str(y)+'_{0:02}'.format(m)+ '_depthi' + str(idepth) + Config['suffixe_all']+'.pickle')<50000:
                have_all_depth = False
    if have_all_depth:
        print('OK :', idepth)