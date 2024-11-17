#!/bin/bash

#SBATCH --nodes=3
#SBATCH --time=24:00:00
#SBATCH --partition=amilan
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --job-name=mtags
#SBATCH --mail-user=cliff.buenodemesquita@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH -C cpu
#SBATCH -q normal
#SBATCH --output=mtags.out

export OMP_NUM_THREADS=32

for i in `ls -1 /scratch/alpine/clbd1748/Australia/QC_reads1/*_R1p.fastq.gz | sed 's/_R1p.fastq.gz//'`
do
mtags profile -f $i\_R1p.fastq.gz -r $i\_R2p.fastq.gz -o /scratch/alpine/clbd1748/Australia/mtags_outputfiles -n $i -t 32 -ma 1000 -mr 1000
done