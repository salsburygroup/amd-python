#!/usr/bin/env python

import numpy as np
import mdtraj
import argparse

parser = argparse.ArgumentParser(
    description='Generate a new trajectory based on selection',
    add_help=False
)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')

inputs.add_argument(
    '-s1',
    action='store',
    dest='structure1',
    help='Structure file corresponding to trajectory',
    type=str,
    required=True
)

inputs.add_argument(
    '-s2',
    action='store',
    dest='structure2',
    help='Structure file corresponding to trajectory',
    type=str,
    required=True
)

inputs.add_argument(
    '-t1',
    action='store',
    dest='trajectory1',
    help='Trajectory1',
    type=str,
    required=True
)

inputs.add_argument(
    '-t2',
    action='store',
    dest='trajectory2',
    help='Trajectory2',
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

inputs.add_argument('-l',
    action='store',
    dest='stride',
    help='stride to use. It will be ignored when frameFile flas is specified.',
    type=int,
    default=1
)

inputs.add_argument(
    '-f',
    action='store',
    dest='frameFile',
    help='file with frame of interest',
    type=str,
    default='noFile'
)

inputs.add_argument(
    '-o',
    action='store',
    dest='outputfile',
    help='output dcd',
    type=str,
    default='stack.dcd'
)

# Parse into useful form
UserInput = parser.parse_args()

if UserInput.frameFile is 'noFile':
    traj1 = mdtraj.load(UserInput.trajectory1, top=UserInput.structure1, stride=UserInput.stride)
    traj2 = mdtraj.load(UserInput.trajectory2, top=UserInput.structure2, stride=UserInput.stride)
else:
    frame = np.loadtxt(UserInput.frameFile,delimiter=',') 
    temp1 = mdtraj.load(UserInput.trajectory1, top=UserInput.structure1)
    traj1 = temp.slice(frame.astype(int)-1, copy = True)
    temp2 = mdtraj.load(UserInput.trajectory2, top=UserInput.structure2)
    traj2 = temp.slice(frame.astype(int)-1, copy = True)

stack = traj1.stack(traj2)
sel = stack.topology.select(UserInput.sel)
extractTraj = stack.atom_slice(sel)

extractTraj.save_dcd(UserInput.outputfile)
