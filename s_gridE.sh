#!/bin/bash
#SBATCH --account=ctb-cdufour
#SBATCH --time=10:0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --array=9

module load python/3.10
module load netcdf
module load mpi4py
#export PYTHONPATH=$PYTHONPATH:
source /home/nplanat/eddytools2/bin/activate

for m in {1..12..1}
do
    python -u /home/nplanat/scripts/Eddies/VF/s_gridE.py $SLURM_ARRAY_TASK_ID ${m}
done