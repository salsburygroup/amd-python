#!/usr/bin/env python
import mdtraj as md
import numpy as np
import argparse
import os

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Calculate, save and plot RMSF', add_help=False)

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
                    default=None)

inputs.add_argument('-l',
                    action='store',
                    dest='stride',
                    help='number of stride',
                    type=str,
                    default=1)

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output prefix for text and png',
                    type=str,
                    required=True)

# Parse into useful form
UserInput = parser.parse_args()

if UserInput.sel==None:
    traj=md.load(UserInput.trajectory, stride=UserInput.stride, top=UserInput.structure)
else:
    traj=md.load(UserInput.trajectory, stride=UserInput.stride, atom_indices=UserInput.sel, top=UserInput.structure)

print(traj.xyz[0])

np.savetxt(UserInput.out_dir, traj.xyz)
