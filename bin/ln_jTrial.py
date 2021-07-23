#! /usr/bin/env python
import argparse
import subprocess

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Compute the transfer entropy for the selected residues', add_help=False)

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
                    default='name CA')

inputs.add_argument('-l',
                    action='store',
                    dest='stride',
                    help='Stride to use',
                    type=str,
                    required=True)

inputs.add_argument('-tau',
                    action='store',
                    dest='tau',
                    help='Decay time (ns)',
                    type=str,
                    default=1)

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output folder for data',
                    type=str,
                    default='transfer_entropy')

# Parse into useful form
UserInput = parser.parse_args()

# Find Helper scripts
script = '/home/wud18/python/ln_j.py'

# Calculating distance
python_cmd = ('sbatch --export=sc=' + script +
              ',s=' + UserInput.structure +
              ',t=' + UserInput.trajectory +
              ',sel=' + UserInput.sel +
              ',l=' + UserInput.stride +
              ',tau=' + UserInput.tau +
              ',o=' + UserInput.out_dir +
              ' /home/wud18/bash/ln_jSubmit.slurm')
print(python_cmd)
subprocess.call(python_cmd, shell=True)