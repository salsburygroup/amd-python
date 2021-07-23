#!/usr/bin/env python
import numpy as np
import mdtraj
import argparse

# Jiajie Xiao
# Jan 14, 2017

parser = argparse.ArgumentParser(
    description='Generate a new trajectory based on selection',
    add_help=False
)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')

inputs.add_argument(
    '-s',
    action='store',
    dest='structure',
    help='Structure file corresponding to trajectory',
    type=str,
    required=True
)

inputs.add_argument(
    '-t',
    action='store',
    dest='trajectory',
    help='Trajectory',
    type=str,
    required=True
)

inputs.add_argument(
    '-sel',
    action='store',
    dest='sel',
    help='selection',
    type=str,
    default='all'
)

inputs.add_argument(
    '-f1',
    action='store',
    dest='N1',
    help='start frame number',
    type=int,
    required=True        
)

inputs.add_argument(
    '-f2',
    action='store',
    dest='N2',
    help='end frame number',
    type=int,
    required=True
)

inputs.add_argument(
    '-o',
    action='store',
    dest='outputfile',
    help='output dcd',
    type=str,
    default='selected.dcd'
)

# Parse into useful form
UserInput = parser.parse_args()
traj = mdtraj.load(UserInput.trajectory, top=UserInput.structure, frame=UserInput.N1-1)

for i in range (UserInput.N1, UserInput.N2):
    temp = mdtraj.load(UserInput.trajectory, top=UserInput.structure, frame=i)
    traj = traj.join(temp)
    
sel = traj.topology.select(UserInput.sel)
extractTraj = traj.atom_slice(sel)

extractTraj.save_dcd(UserInput.outputfile)
