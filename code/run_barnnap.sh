#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --partition=amilan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=barnnap
#SBATCH --mail-user=cliff.buenodemesquita@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH -C cpu
#SBATCH -q normal
#SBATCH --output=barnnap.out

export OMP_NUM_THREADS=16
module purge
module load python
source /home/clbd1748/.bashrc
source activate barnnap_env

# Process all genomes to extract 16S rRNA sequences using Barrnap
for genome in /scratch/alpine/clbd1748/Australia/BradyStrains/*.fna; do

    # Run Barrnap
    barrnap "$genome" --kingdom 'bac' --quiet --threads 16 --outseq "${genome%.fasta}_16S.fasta" > "${genome%.fasta}_16S.gff"
    
done