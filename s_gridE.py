import sys
import f_gridE



y = int(eval(sys.argv[1]))
m = int(eval(sys.argv[2]))


from Config_file_REF12_3 import Config

f_gridE.get_gridE(y,m,Config)