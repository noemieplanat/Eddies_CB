#!/bin/bash
#SBATCH --account=ctb-cdufour
#SBATCH --time=15:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=5G
#SBATCH --array=9

module load python/3.10
module load netcdf
module load mpi4py
#export PYTHONPATH=$PYTHONPATH:
source /home/nplanat/eddytools2/bin/activate

python -u /home/nplanat/scripts/Eddies/VF_3/s_sensi.py  $SLURM_ARRAY_TASK_ID
