#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --partition=amilan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=align
#SBATCH --mail-user=cliff.buenodemesquita@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH -C cpu
#SBATCH -q normal
#SBATCH --output=align.out

export OMP_NUM_THREADS=16
module purge
module load python
source /home/clbd1748/.bashrc
source activate gtdbtk_env

gtdbtk align --identify_dir /scratch/alpine/clbd1748/Australia/identify_output/ --out_dir /scratch/alpine/clbd1748/Australia/align_output --cpus 16 --skip_gtdb_refs --tmpdir /scratch/alpine/clbd1748/Australia/gtdb_tmp