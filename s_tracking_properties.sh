#!/bin/bash
#SBATCH --account=ctb-cdufour
#SBATCH --time=3:0:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --array=0-49

module load python/3.10
module load netcdf
module load mpi4py
#export PYTHONPATH=$PYTHONPATH:
source /home/nplanat/eddytools2/bin/activate

python -u /home/nplanat/scripts/Eddies/VF_3/s_tracking_properties.py $SLURM_ARRAY_TASK_ID