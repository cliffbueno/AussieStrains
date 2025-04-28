#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=6:00:00
#SBATCH --partition=amilan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=raxml
#SBATCH --mail-user=cliff.buenodemesquita@colorado.edu
#SBATCH --mail-type=ALL
#SBATCH -C cpu
#SBATCH -q normal
#SBATCH --output=raxml.out

export OMP_NUM_THREADS=16
module purge
module load python
source /home/clbd1748/.bashrc
source activate raxml_env

raxmlHPC -f a -m PROTGAMMALG -p 12345 -x 12345 -N 1000 -s alignment.phy -n LG_BOOTSTRAP -o Nitrobacter.fna