#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition=amilan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=pogenom
#SBATCH --mail-user=cliff.buenodemesquita@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH -C cpu
#SBATCH -q normal
#SBATCH --output=pogenom.out

export OMP_NUM_THREADS=16
module purge
module load python
source /home/clbd1748/.bashrc

# Run POGENOM
perl pogenom.pl --vcf_file /scratch/alpine/clbd1748/POGENOM/Input_POGENOM/06_VCF/Aus/params_cov_0_bdth_0_subsamp_FALSE_mpq_20_bq_15/BradyReference.vcf --out Brady --gff_file BradyReference_fixed.gff --fasta_file BradyReference.fasta --genetic_code_file standard_genetic_code.txt