#!/bin/bash
#SBATCH --account=ctb-cdufour
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --array=14,20,27,32,39

module load python/3.10
module load netcdf
module load mpi4py
#export PYTHONPATH=$PYTHONPATH:
source /home/nplanat/eddytools2/bin/activate

python -u /home/nplanat/scripts/Eddies/VF_3/s_format.py $SLURM_ARRAY_TASK_ID 0
