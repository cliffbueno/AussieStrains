#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --partition=amilan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name=coverm_brady
#SBATCH --mail-user=cliff.buenodemesquita@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH -C cpu
#SBATCH -q normal
#SBATCH --output=coverm_brady.out

export OMP_NUM_THREADS=32

for i in `ls -1 /scratch/alpine/clbd1748/Australia/QC_reads/*_R1p.fastq.gz | sed 's/_R1p.fastq.gz//'`
do
coverm genome -c $i\_R1p.fastq.gz $i\_R2p.fastq.gz -d /scratch/alpine/clbd1748/Australia/brady -x fna -m relative_abundance -o output.tsv --output-format dense -t 32
done


