#! /usr/bin/env python

#####
#Principal_axes.py   
#Mon Jan 25 11:56:47 EST 2016
#Ryan Melvin
#####

import MDAnalysis as md
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import seaborn


# Initialize parser
parser = argparse.ArgumentParser(description = 'Calculate Principal axes of selection', add_help=False) 

#List all possible user input
inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-top', action='store', dest='top',help='topology file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-traj', action='store', dest='traj',help='Trajectory',type=str,required=True)
inputs.add_argument('-sel', action='store', dest='sel', help='Atom selection',type=str,default='all')
inputs.add_argument('-o', action='store', dest='out_name',help='Output files prefix ',type=str,required=True)

UserInput = parser.parse_args()


# Load Trajectory
u = md.Universe(UserInput.top,UserInput.traj)

# Select atoms
group = u.select_atoms(UserInput.sel)

# Setup storage data frames
axes = pd.DataFrame(columns = ['e1x', 'e1y', 'e1z','e2x', 'e2y', 'e2z','e3x', 'e3y', 'e3z'])
dots = pd.DataFrame(columns = ['e1','e2','e3'])

# Get principal axes of first frame for reference
group.align_principal_axis(0,[0.,1.,0.])
ref = group.principal_axes()

for ts in u.trajectory:
    group.align_principal_axis(0,[0.,1.,0.])
    e = group.principal_axes()
    if np.round(np.dot(ref[0],np.transpose(e[0]))) == -1:
        e[0] = -e[0]
        e[1] = -e[1]
    axes.loc[ts.frame] = np.concatenate((e[0],e[1],e[2]), axis = 0)
    # Calculate dot products
    dots.loc[ts.frame] = [
            np.dot(ref[0],np.transpose(e[0])),
            np.dot(ref[1],np.transpose(e[1])),
            np.dot(ref[2],np.transpose(e[2])),
            ]

axes.to_csv(UserInput.out_name + '_axes.txt', sep='\t', index=False)
dots.to_csv(UserInput.out_name + '_dots.txt', sep='\t', index=False)

t = np.arange(len(u.trajectory)) #make a list of frames

e1plot = plt.scatter(t, dots['e1'], color='red')
e2plot = plt.scatter(t, dots['e2'], color='blue')
e3plot = plt.scatter(t, dots['e3'], color='green')

plt.legend((e1plot, e2plot, e3plot),
        ('e1', 'e2', 'e3'),
        scatterpoints=1,
        loc='upper right',
        ncol=3,
        fontsize=10)
plt.savefig(UserInput.out_name + '_dots.png')
