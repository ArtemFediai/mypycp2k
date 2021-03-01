#!/bin/bash 
#
#SBATCH --time=00:35:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --partition=batch
#SBATCH --job-name=interessante_01

# first non-empty non-comment line ends SLURM options

module load cp2k/7.1/

mpirun -np 20 $CP2K -i DIAG_GW_PBE_f6tcnnq.inp -o out.out

date
