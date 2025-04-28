#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition=amilan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name=variantcalling
#SBATCH --mail-user=cliff.buenodemesquita@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH -C cpu
#SBATCH -q normal
#SBATCH --output=variantcalling.out

export OMP_NUM_THREADS=32
module purge
module load python
source /home/clbd1748/.bashrc
source activate variant_env

# Run the variant calling pipeline
/scratch/alpine/clbd1748/Australia/variant_calling_pipeline.sh