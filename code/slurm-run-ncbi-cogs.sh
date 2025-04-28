#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --partition=amilan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=anvio-run-ncbi-cogs
#SBATCH --mail-user=cliff.buenodemesquita@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH -C cpu
#SBATCH -q normal
#SBATCH --output=anvio-run-ncbi-cogs.out

export OMP_NUM_THREADS=16
module purge
module load python
source /home/clbd1748/.bashrc
source activate anvio-8

# Loop through each -contig database and annotate with COGs
for db in *.db; do anvi-run-ncbi-cogs -c "$db" --num-threads 16; done