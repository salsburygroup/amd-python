#! /usr/bin/env python
import argparse
import subprocess

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='PALD clustering method', add_help=False)

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
                    default='protein and resid 167:170')

inputs.add_argument('-title',
                    action='store',
                    dest='title',
                    help='Title of the plot',
                    type=str,
                    default='Closest mean distance between Na$\mathregular{^+}$ and 220s loop')

inputs.add_argument('-tm',
                    action='store',
                    dest='timestep',
                    help='time interval(ps)',
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
script = '/home/wud18/python/pald.py'

# Calculating distance
python_cmd = ('sbatch --export=sc=' + script +
              ',s=' + UserInput.structure +
              ',t=' + UserInput.trajectory +
              ',sel=' + UserInput.sel +
              ',title=' + UserInput.title +
              ',tm=' + UserInput.timestep +
              ',o=' + UserInput.out_dir +
              ' /home/wud18/bash/paldSubmit.slurm')
print(python_cmd)
subprocess.call(python_cmd, shell=True)