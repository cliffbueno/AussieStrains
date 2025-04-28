#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --partition=amilan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=anvio-contigs-db
#SBATCH --mail-user=cliff.buenodemesquita@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH -C cpu
#SBATCH -q normal
#SBATCH --output=anvio-contigs-db.out

export OMP_NUM_THREADS=16
module purge
module load python
source /home/clbd1748/.bashrc
source activate anvio-8

# Loop through each -contigs-fasta artifact and make a contigs-db
for file in *-contigs-fasta; do anvi-gen-contigs-database -f "$file" --project-name "${file%-contigs-fasta}" -o "${file%-contigs-fasta}.db"; done