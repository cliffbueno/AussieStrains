#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --partition=amilan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=sylph_profile_genomes_ani95.sh
#SBATCH --mail-user=cliff.buenodemesquita@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH -C cpu
#SBATCH -q normal
#SBATCH --output=sylph_profile_genomes_ani95.out

export OMP_NUM_THREADS=16

module purge
module load python
source /home/clbd1748/.bashrc
source activate sylph_env

sylph profile /scratch/alpine/clbd1748/Australia/ref_genomes_sketch.syldb /scratch/alpine/clbd1748/Australia/reads_sketch/*.sylsp -t 16 -m 95 > /scratch/alpine/clbd1748/Australia/sylph_profile_ani95.tsv