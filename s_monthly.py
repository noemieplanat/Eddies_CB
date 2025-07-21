import sys
import f_monthly
import xarray as xr
import numpy as np


depth_index = int(eval(sys.argv[1]))
from Config_file_REF12_3 import Config, Model  

yi = 1995 ; yf = 2020 
f_monthly.make_ORCA_monthly_m_filter_dt(yi, yf, depth_index, Config, filter_time = 2)
f_monthly.make_ORCA_monthly_m_filter_prop(yi, yf, depth_index, Config, filter_time = 2)