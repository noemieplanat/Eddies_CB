#!/bin/bash
#SBATCH --account=ctb-cdufour
#SBATCH --time=10:00
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=1

module load netcdf
module load nco


REF12_a='/home/nplanat/projects/ctb-cdufour/shared_files/models/CREG12-REF12/1a'
REF12_m='/home/nplanat/projects/ctb-cdufour/shared_files/models/CREG12-REF12/1m'
PVC='/home/nplanat/projects/ctb-cdufour/shared_files/models/CREG12-REF12/climato'


#ncra -O -4 -v mldr1_1 $REF12_a/CREG12.L75-REF12_y*.1d_gridT.nc $PVC/MLD_1995_2020.nc
#ncra -O -4 -v mldr1_1 $REF12_m/CREG12.L75-REF12_y*m0[1-3].1d_gridT.nc $PVC/MLD_JFM_1995_2020.nc
#ncra -O -4 -v mldr1_1 $REF12_m/CREG12.L75-REF12_y*m0[4-6].1d_gridT.nc $PVC/MLD_AMJ_1995_2020.nc
#ncra -O -4 -v mldr1_1 $REF12_m/CREG12.L75-REF12_y*m0[7-9].1d_gridT.nc $PVC/MLD_JAS_1995_2020.nc
ncra -O -4 -v mldr1_1 $REF12_m/CREG12.L75-REF12_y*m1[0-2].1d_gridT.nc $PVC/MLD_OND_1995_2020.nc