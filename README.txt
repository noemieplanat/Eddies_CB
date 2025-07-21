All steps cleaned for detection and tracking with constant std_OW 
REF12 version 3

1. Interpolation [per depth index and per year]
    Compute OW, vorticity and std_OW from daily outputs and interp t/f points
    * s_interp.sh -> s_interp.py -> f_interp.py 
    average all OW_std per depth level into one file [bridge REF12_2 onto REF12_3]
    * s_climato_OW.sh YTD 
2. detection
    Use the mean std from s_climato_OW to detect eddies. 
    * s_detect.sh -> s_detect.py -> f_detect.py

3. tracking : I track also eddies that live only for one day (separately)
    * s_tracking.sh -> s_tracking.py -> f_tracking.py
4. format : 
    Format into netcdf filtered to the region and duration 
    *s_format.sh -> s_format.py -> f_tracking.py
    Format into 2D fields
    *s_monthly.sh -> s_monthly.py -> f_monthly.py

    