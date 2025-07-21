#!/bin/bash
#SBATCH --account=ctb-cdufour
#SBATCH --time=8:0:0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --array=1995

module load python/3.10
module load netcdf
module load mpi4py
#export PYTHONPATH=$PYTHONPATH:
source /home/nplanat/eddytools2/bin/activate
for i in {0..0..1}
do
    echo "Starting depth  ${i}"
    python /home/nplanat/scripts/Eddies/VF_3/s_interp.py ${i} $SLURM_ARRAY_TASK_ID Config_file_G12
done
