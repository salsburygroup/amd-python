# coding: utf-8
#!/usr/bin/env python
import numpy as np
import MDAnalysis
import argparse

# Based on material at https://groups.google.com/forum/#!topic/mdnalysis-discussion/ub6tBmZGpbw

# Jiajie Xiao
# Oct 18, 2016

parser = argparse.ArgumentParser(
    description='Output the mean structure in a trajectory',
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

inputs.add_argument(
    '-outputfile',
    action='store',
    dest='outputfile',
    help='output mean structure pdb',
    type=str,
    default='meanStructure.pdb'
)


# Parse into useful form
UserInput = parser.parse_args()

u = MDAnalysis.Universe(UserInput.structure, UserInput.trajectory)
p = u.select_atoms(UserInput.sel)

p_avg = np.zeros_like(p.positions)

for i in u.trajectory:
	p_avg += p.positions
	
p_avg /= len(u.trajectory)
    
p.set_positions(p_avg)
#pdb = MDAnalysis.Writer(UserInput.outputfile, multiframe=False)
#pdb.write(p)
#pdb.close()
p.write(UserInput.outputfile)	
# in my test, the output pdb has a binary comment showing which plugin was used for writing the pdb