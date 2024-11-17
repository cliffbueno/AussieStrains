#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --partition=amilan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=sylph_sketch_genomes
#SBATCH --mail-user=cliff.buenodemesquita@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH -C cpu
#SBATCH -q normal
#SBATCH --output=sylph_sketch_genomes.out

export OMP_NUM_THREADS=16

module purge
module load python
source /home/clbd1748/.bashrc
source activate sylph_env

sylph sketch -g /scratch/alpine/clbd1748/Australia/Ref_genomes/*.fna.gz -o /scratch/alpine/clbd1748/Australia/ref_genomes_sketch -t 16