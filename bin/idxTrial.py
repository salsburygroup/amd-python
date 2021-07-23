#! /usr/bin/env python
import argparse
import subprocess
import os

print('Make sure to use the absolute path for output directory(-o)')
# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Compute the vector difference between Na+ binding loop and Na+', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-s',
                    action='store',
                    dest='structure',
                    help='Structure file corresponding to trajectory',
                    type=str,
                    required=True)

inputs.add_argument('-t',
                    action='store',
                    dest='trajectory',
                    help='Trajectory',
                    type=str,
                    required=True)

inputs.add_argument('-sel',
                    action='store',
                    dest='sel',
                    help='Atom selection',
                    type=str,
                    default='protein\ and\ resid\ 264:274')

inputs.add_argument('-dist',
                    action='store',
                    dest='dist',
                    help='Threshold for the cutoff distance',
                    type=str,
                    required=True)

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output folder for data',
                    type=str,
                    required=True)

# Parse into useful form
UserInput = parser.parse_args()

# Find Helper scripts
script = '/home/wud18/python/idx.py'

# Calculating distance
if not os.path.exists(UserInput.out_dir):
    os.makedirs(UserInput.out_dir)
python_cmd = ('sbatch --export=sc=' + script +
              ',s=' + UserInput.structure +
              ',t=' + UserInput.trajectory +    
              ',sel=' + UserInput.sel + 
              ',dist=' + UserInput.dist + 
              ',o=' + UserInput.out_dir + 
              ' /home/wud18/bash/idxSubmit.slurm')
print(python_cmd)
subprocess.call(python_cmd, shell=True)
