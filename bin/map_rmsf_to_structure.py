import pandas as pd
import numpy as np
import os
import argparse

# Initialize parser.
parser = argparse.ArgumentParser(description='Map rmsf into structure ', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument(
    '-s',
    action='store',
    dest='structure',
    help='Structure file corresponding to trajectory',
    type=str,
    default='protein.pdb'
)

inputs.add_argument(
    '-d',
    action='store',
    dest='d',
    help='Rmsf of the protein',
    type=str,
    default='rmsf_averaged_protein/rmsf_averaged_protein.dat'
)

inputs.add_argument(
    '-o',
    action='store',
    dest='out_dir',
    help='Output directory',
    type=str,
    default='rmsf_map.pdb'
)

# Parse into useful form
UserInput = parser.parse_args()

# Make output directory
fname = './rmsf2structure/'
if not os.path.exists(fname):
    os.mkdir(fname)

# open data
path=os.getcwd()
data = open(path + '/' + UserInput.d)
structure = open(path + '/' + UserInput.structure)

# read data line by line
data_new = [1,1,1]
for line in data.readlines():
    data_new.append(float(line.strip()))
data_new = data_new + [1,1,1]

# replace occupency column with rmsf
n = 0
f=open(fname + UserInput.out_dir, 'w')
for line in structure:
    print(line.strip().replace('  0.00  ', '  ' + format(data_new[n],'.4f') + '  '), file=f)
    n = n+1
