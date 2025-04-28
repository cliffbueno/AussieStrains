#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=6:00:00
#SBATCH --partition=amilan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=fastani
#SBATCH --mail-user=cliff.buenodemesquita@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH -C cpu
#SBATCH -q normal
#SBATCH --output=fastani.out

export OMP_NUM_THREADS=16
module purge
module load python
source /home/clbd1748/.bashrc
source activate fastani_env

fastANI --ql /scratch/alpine/clbd1748/Australia/queryList.txt --rl /scratch/alpine/clbd1748/Australia/referenceList.txt -o /scratch/alpine/clbd1748/Australia/Brady_fastani.txt