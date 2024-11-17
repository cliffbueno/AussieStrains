#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --partition=amilan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=sylph_query_genomes_ani90
#SBATCH --mail-user=cliff.buenodemesquita@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH -C cpu
#SBATCH -q normal
#SBATCH --output=sylph_query_genomes_ani90.out

export OMP_NUM_THREADS=16

module purge
module load python
source /home/clbd1748/.bashrc
source activate sylph_env

sylph query /scratch/alpine/clbd1748/Australia/ref_genomes_sketch.syldb /scratch/alpine/clbd1748/Australia/reads_sketch/*.sylsp -t 16 -m 90 > /scratch/alpine/clbd1748/Australia/sylph_query_ani90.tsv