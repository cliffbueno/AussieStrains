#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition=amilan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=anvio-pan-genome
#SBATCH --mail-user=cliff.buenodemesquita@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH -C cpu
#SBATCH -q normal
#SBATCH --output=anvio-pan-genome.out

export OMP_NUM_THREADS=16
module purge
module load python
source /home/clbd1748/.bashrc
source activate anvio-8

# Loop through each -contig database and annotate with KOs
anvi-pan-genome -g BRADY-GENOMES.db \
                --project-name "Bradyrhizobium_Pan" \
                --output-dir BRADY \
                --num-threads 16 \
                --minbit 0.5 \
                --mcl-inflation 10 \
                --use-ncbi-blast