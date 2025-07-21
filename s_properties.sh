#!/bin/bash
#SBATCH --account=ctb-cdufour
#SBATCH --time=5:00:0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30G
#SBATCH --array=1995-2020

module load python/3.10
module load netcdf
module load mpi4py
#export PYTHONPATH=$PYTHONPATH:
source /home/nplanat/eddytools2/bin/activate

for de in {0..49..1}
do
    for m in {1..12..1}
    do
        python /home/nplanat/scripts/Eddies/VF_3/s_properties.py $SLURM_ARRAY_TASK_ID ${de} ${m}
    done
done
