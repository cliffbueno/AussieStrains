#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --partition=amilan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=sylph_sketch_reads
#SBATCH --mail-user=cliff.buenodemesquita@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH -C cpu
#SBATCH -q normal
#SBATCH --output=sylph_sketch_reads.out

export OMP_NUM_THREADS=16

module purge
module load python
source /home/clbd1748/.bashrc
source activate sylph_env

sylph sketch -1 /scratch/alpine/clbd1748/Australia/QC_reads/*R1p.fastq.gz -2 /scratch/alpine/clbd1748/Australia/QC_reads/*R2p.fastq.gz -d /scratch/alpine/clbd1748/Australia/reads_sketch -t 16