#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --partition=amilan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name=coverm_top4
#SBATCH --mail-user=cliff.buenodemesquita@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH -C cpu
#SBATCH -q normal
#SBATCH --output=coverm_top4.out

export OMP_NUM_THREADS=32

module purge
module load python
source /home/clbd1748/.bashrc
source activate coverm_env

for i in `ls -1 /scratch/alpine/clbd1748/Australia/QC_reads/*_R1p.fastq.gz | sed 's/_R1p.fastq.gz//'`
do
coverm genome -c $i\_R1p.fastq.gz $i\_R2p.fastq.gz -d /scratch/alpine/clbd1748/Australia/top4 -x fna.gz -m relative_abundance -o output_top4.tsv --output-format dense -t 32
done


