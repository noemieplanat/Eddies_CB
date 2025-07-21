#!/bin/bash
#SBATCH --account=ctb-cdufour
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=1

module load netcdf
module load nco


PVI='/home/nplanat/scratch/REF12/Interp/Interp_CREG_1d_std50_' #bridge version REF12_2 onto REF12_3
PVM='/home/nplanat/projects/ctb-cdufour/shared_files/models/CREG12-REF12/1m'
PVC='/home/nplanat/projects/ctb-cdufour/shared_files/models/CREG12-REF12/climato'

#change time to record dimension
for d in {0..13..1}
do
    for y in {1995..2020..1}
    do
        for m in {1..12..1}
        do
            case $m in
                [0-9])
                    m="0$m"
                    ;;
                *)
            esac
            ncks -v OW_std $PVI${y}_${m}_depthi${d}REF12_0.nc $PVM/OW_std_${y}_${m}_depthi${d}.1m.nc
            ncecat -O -u time $PVM/OW_std_${y}_${m}_depthi${d}.1m.nc $PVM/OW_std_${y}_${m}_depthi${d}.1m.nc
        done
    done
    ncra -O -4 -v OW_std $PVM/OW_std_*_depthi${d}.1m.nc $PVC/OW_depthi${d}.1995_2020.nc
done

