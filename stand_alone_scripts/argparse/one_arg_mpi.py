#!/home/artem/anaconda3/bin/python
##!/bin/bash/ python3
import argparse
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

parser = argparse.ArgumentParser(description='Return the argument')
parser.add_argument('my_argument',
                    help='my_argument')

args = parser.parse_args()


if rank == 0:
    data = {'a': args.my_argument}
    comm.send(data, dest=1, tag=11)
    print(f"this is from rank 0. I try to send the following: {data}")
elif rank == 1:
    data = comm.recv(source=0, tag=11)
    print(f"this is from rank 1. I have recieved the following: {data}")

print(f'the following argument is passed: {args.my_argument}')
