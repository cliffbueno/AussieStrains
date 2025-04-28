#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --partition=amilan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=prottest
#SBATCH --mail-user=cliff.buenodemesquita@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH -C cpu
#SBATCH -q normal
#SBATCH --output=prottest.out

export OMP_NUM_THREADS=16
module purge
module load python
source /home/clbd1748/.bashrc

java -jar /scratch/alpine/clbd1748/Australia/prottest-3.4.2/prottest-3.4.2.jar -i /scratch/alpine/clbd1748/Australia/align_output/align/gtdbtk.bac120.user_msa.fasta -all-matrices -all-distributions -threads 16