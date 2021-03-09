#!/home/artem/anaconda3/bin/python
##!/bin/bash/ python3
import argparse

parser = argparse.ArgumentParser(description='Return the argument')
parser.add_argument('my_argument',
                    help='my_argument')

args = parser.parse_args()


print(f'the following argument is passed: {args.my_argument}')
