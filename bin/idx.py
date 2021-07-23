#! /usr/bin/env python
import os
import MDAnalysis
import argparse
import numpy as np
from scipy.spatial import distance

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
                    default='protein and resid 264:274')

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

u = MDAnalysis.Universe(os.path.join(UserInput.out_dir, UserInput.structure), os.path.join(UserInput.out_dir, UserInput.trajectory))
sodiumLoop=u.select_atoms(UserInput.sel)
sodiumGrp=u.select_atoms('resname SOD')
T=len(u.trajectory)
data=np.zeros((T, len(sodiumLoop)*3))
data_no_noise=np.zeros((T, len(sodiumLoop)*3))
index=[]
for i in range(T):
    u.trajectory[i] #Go to current frame
    Na_dist_mean=np.mean(distance.cdist(sodiumLoop.positions,sodiumGrp.positions,'euclidean'),axis=0)
    Na_min_index=np.argmin(Na_dist_mean)
    a=sodiumGrp.positions[Na_min_index]
    b=sodiumLoop.positions
    data[i]=(a-b).flatten()    
    if np.min(Na_dist_mean)>float(UserInput.dist):
        data_no_noise[i]=(a-b).flatten()
        index=index + [i+1]
np.savetxt(os.path.join(UserInput.out_dir, 'idx' + UserInput.dist), index, delimiter=' ')
