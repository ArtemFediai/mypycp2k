#!/bin/bash 
#
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=20
#SBATCH --partition=batch
#SBATCH --job-name=interesante

# first non-empty non-comment line ends SLURM options

module load cp2k/7.1/
#module load gnu8/8.3.0 openmpi3
#module load devel/python/3.6
#module load mpi/impi/2020

$conda activate cp2k_environment

python3 extract_homo_lumo_and_gw_energies.py > out.out
#srun python3 extract_homo_lumo_and_gw_energies.py > out.out

date
