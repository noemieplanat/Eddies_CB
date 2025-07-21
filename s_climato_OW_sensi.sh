#!/bin/bash
#SBATCH --account=ctb-cdufour
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=1

module load netcdf
module load nco


PVI='/home/nplanat/scratch/REF12/Interp/' #bridge version REF12_2 onto REF12_3
PVI3='/home/nplanat/scratch/REF12_3/Interp/' 
PVM='/home/nplanat/projects/ctb-cdufour/shared_files/models/CREG12-REF12/1m'
PVC='/home/nplanat/projects/ctb-cdufour/shared_files/models/CREG12-REF12/climato'


REF12_3/Interp/Interp_CREG_1d_std100_2000_01_depthi14REF12_3.nc
#change time to record dimension
for d in {14,20,27,32,29}
do
    y=2000
    for m in {1..12..1}
    do
        case $m in
            [0-9])
                m="0$m"
                ;;
            *)
        esac
        #ncks -v OW_std ${PVI}Interp_CREG_1d_std50_${y}_${m}_depthi${d}REF12_0.nc $PVM/OW_2000_std_${y}_${m}_depthi${d}.1m.nc
        ncks -v OW_std ${PVI3}Interp_CREG_1d_std12_${y}_${m}_depthi${d}REF12_3.nc $PVM/OW_sig12_std_${y}_${m}_depthi${d}.1m.nc
        ncks -v OW_std ${PVI3}Interp_CREG_1d_std24_${y}_${m}_depthi${d}REF12_3.nc $PVM/OW_sig24_std_${y}_${m}_depthi${d}.1m.nc
        ncks -v OW_std ${PVI3}Interp_CREG_1d_std74_${y}_${m}_depthi${d}REF12_3.nc $PVM/OW_sig74_std_${y}_${m}_depthi${d}.1m.nc
        ncks -v OW_std ${PVI3}Interp_CREG_1d_std100_${y}_${m}_depthi${d}REF12_3.nc $PVM/OW_sig100_std_${y}_${m}_depthi${d}.1m.nc
        #ncecat -O -u time $PVM/OW_std_${y}_${m}_depthi${d}.1m.nc $PVM/OW_std_${y}_${m}_depthi${d}.1m.nc
        ncecat -O -u time $PVM/OW_sig12_std_${y}_${m}_depthi${d}.1m.nc $PVM/OW_sig12_std_${y}_${m}_depthi${d}.1m.nc
        ncecat -O -u time $PVM/OW_sig24_std_${y}_${m}_depthi${d}.1m.nc $PVM/OW_sig24_std_${y}_${m}_depthi${d}.1m.nc
        ncecat -O -u time $PVM/OW_sig74_std_${y}_${m}_depthi${d}.1m.nc $PVM/OW_sig74_std_${y}_${m}_depthi${d}.1m.nc
        ncecat -O -u time $PVM/OW_sig100_std_${y}_${m}_depthi${d}.1m.nc $PVM/OW_sig100_std_${y}_${m}_depthi${d}.1m.nc

    done
    ncra -O -4 -v OW_std $PVM/OW_std_2000*_depthi${d}.1m.nc $PVC/OW_depthi${d}.2000.nc
    ncra -O -4 -v OW_std $PVM/OW_sig12_std_2000*_depthi${d}.1m.nc $PVC/OW_sig12_std_${d}.2000.nc
    ncra -O -4 -v OW_std $PVM/OW_sig24_std_2000*_depthi${d}.1m.nc $PVC/OW_sig24_std_${d}.2000.nc
    ncra -O -4 -v OW_std $PVM/OW_sig74_std_2000*_depthi${d}.1m.nc $PVC/OW_sig74_std_${d}.2000.nc
    ncra -O -4 -v OW_std $PVM/OW_sig100_std_2000*_depthi${d}.1m.nc $PVC/OW_sig100_std_${d}.2000.nc
done

