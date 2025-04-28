#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition=amilan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --job-name=annotate-contigs
#SBATCH --mail-user=cliff.buenodemesquita@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH -C cpu
#SBATCH -q normal
#SBATCH --output=annotate-contigs.out

export OMP_NUM_THREADS=48
module purge
module load python
source /home/clbd1748/.bashrc
source activate anvio-8

# Loop through each -contig database and run hmms and scan trnas
for db in *.db; do anvi-run-hmms -c "$db" --num-threads 48 --also-scan-trnas; done

# Loop through each -contig database and run scg taxonomy
for db in *.db; do anvi-run-scg-taxonomy -c "$db" --num-threads 48; done

# Loop through each -contig database and annotate with COGs
for db in *.db; do anvi-run-ncbi-cogs -c "$db" --num-threads 48; done

# Loop through each -contig database and annotate with KOs
for db in *.db; do anvi-run-kegg-kofams -c "$db" --num-threads 48; done

# Loop through each -contig database and annotate with CAZymes
for db in *.db; do anvi-run-cazymes -c "$db" --num-threads 48; done