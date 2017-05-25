#!/usr/bin/env python
# coding: utf-8
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
    '-structure',
    action='store',
    dest='structure',
    help='Structure file corresponding to trajectory',
    type=str,
    required=True
)

inputs.add_argument(
    '-traj',
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
inputs.add_argument('-stride',
                    action='store',
                    dest='stride',
                    help='stride to use. It will be ignored when frameFile flas is specified.',
                    type=int,
                    default=1
                    )                   
inputs.add_argument(
    '-frameFile',
    action='store',
    dest='frameFile',
    help='file with frame of interest',
    type=str,
    default='noFile'
)
inputs.add_argument(
    '-outputfile',
    action='store',
    dest='outputfile',
    help='output dcd',
    type=str,
    default='selected.dcd'
)


# Parse into useful form
UserInput = parser.parse_args()

if UserInput.frameFile is 'noFile':
    traj = mdtraj.load(UserInput.trajectory, top=UserInput.structure, stride=UserInput.stride)
else:
    temp = mdtraj.load(UserInput.trajectory, top=UserInput.structure)
    frame = np.loadtxt(UserInput.frameFile,delimiter=',') 
    traj = temp.slice(frame.astype(int)-1, copy = True)
    
sel = traj.topology.select(UserInput.sel)
extractTraj = traj.atom_slice(sel)

extractTraj.save_dcd(UserInput.outputfile)
