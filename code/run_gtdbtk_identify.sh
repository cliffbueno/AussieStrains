#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --partition=amilan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=identify
#SBATCH --mail-user=cliff.buenodemesquita@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH -C cpu
#SBATCH -q normal
#SBATCH --output=identify.out

export OMP_NUM_THREADS=16
module purge
module load python
source /home/clbd1748/.bashrc
source activate gtdbtk_env

gtdbtk identify --genome_dir /scratch/alpine/clbd1748/Australia/BradyStrainsGTDB/ --out_dir /scratch/alpine/clbd1748/Australia/identify_output --cpus 16 -x gz --temp_dir /scratch/alpine/clbd1748/Australia/gtdbtk_temp