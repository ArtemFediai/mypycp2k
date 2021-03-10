#!/home/artem/anaconda3/bin/python

import argparse
import os
from mpi4py import MPI
from cp2k_run.cp2k_run import cp2k_run

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


parser = argparse.ArgumentParser(description='rank and num of cpus')
parser.add_argument('-rank')
parser.add_argument('-num_cpus')

args = parser.parse_args()

print(f'the following argument is passed: {args.rank}')

print('Now I will be doing cp2k...')

run_folder = f'{args.rank}'  # will be created to copy there input file and run cp2k

if not os.path.exists(run_folder):
    os.mkdir(run_folder)

cp2k_run(cp2k_executable='cp2k.popt',
         run_type='mpi',
         np=args.num_cpus,
         xyz_file='input.inp',
         output_file=f'out_{args.rank}.out',
         error_file=f'err_{args.rank}.err',
         execution_directory=f'{args.rank}',
         )

