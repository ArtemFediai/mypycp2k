#!/bin/bash 
#SBATCH --time=02:00:00
#SBATCH --output=./out/out_%A_%a.out
#SBATCH --error=./out/out_%A_%a.err
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=2000
#SBATCH --tasks-per-node=16
#SBATCH --partition=batch
#SBATCH --job-name=104th
#SBATCH --array=1-2000%50
#SBATCH --exclude=bionano01,bionano02,bionano03,bionano04,bionano05,bionano06

PLUS=103000
NEW_SLURM_ARRAY_TASK_ID=$(( $SLURM_ARRAY_TASK_ID + $PLUS ))
echo $NEW_SLURM_ARRAY_TASK_ID

# first non-empty non-comment line ends SLURM options
# cp2k.8 is in your PATH

module load gnu8/8.3.0
module load openmpi3

CP2K_popt=cp2k.popt


cd ${SLURM_SUBMIT_DIR}
echo slurm_submit_dir, $SLURM_SUBMIT_DIR
export SCRATCH=/scratch

mkdir $SCRATCH/bh5670  # if no such folder --> scratch
mkdir $SCRATCH/bh5670/sim  # simulations --> scratch

export OMP_NUM_THREADS=1
export OMP_PROC_BIND=FALSE
export SLURM_CPU_BIND=none
unset SLURM_CPU_BIND

python ~/mypycp2k/examples/234_RI_PBE_yaml_factory.py -rank $NEW_SLURM_ARRAY_TASK_ID -num_cpus $SLURM_NTASKS -i input_factory.yaml
