#! /usr/bin/env python
from Analysis import Distance, Plotter, Saver, TrajectoryReader, TrajectoryProcessor
from scipy.spatial import distance
import matplotlib.pyplot as plt
import mdtraj as md
import pandas as pd
import numpy as np
import argparse
import os

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
    description='Compute (min,max,mean or all pairwise) distance of selected pairs of atoms for given MD trajectories', add_help=False)

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

inputs.add_argument('-l',
                    action='store',
                    dest='stride',
                    help='stride to use',
                    type=int,
                    default=1)

inputs.add_argument('-periodic',
                    action='store',
                    dest='periodic',
                    help='periodic true of false',
                    type=bool,
                    default=True)

inputs.add_argument('-sel1',
                    action='store',
                    dest='sel1',
                    help='selection 1',
                    type=str,
                    required=True)

inputs.add_argument('-sel2',
                    action='store',
                    dest='sel2',
                    help='selection 2',
                    type=str,
                    required=True)

inputs.add_argument('-idx',
                    action='store',
                    dest='idx',
                    help='Output index',
                    type=str,
                    required=True)

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output folder for data',
                    type=str,
                    required=True)

inputs.add_argument('-title',
                    action='store',
                    dest='title',
                    help='Title of the plot',
                    type=str,
                    default='Mean distance between two dimers')

# Parse into useful form
UserInput = parser.parse_args()
traj = md.load(UserInput.trajectory, top=UserInput.structure, stride=UserInput.stride)
atom_pairs = traj.topology.select_pairs(selection1=UserInput.sel1,selection2=UserInput.sel2)

# Make output directory
out_dir=UserInput.out_dir
fname = './' + 'distance_' + out_dir
if not os.path.exists(fname):
    os.mkdir(fname)

# Compute distance
d = md.compute_distances(traj,atom_pairs,periodic=UserInput.periodic,opt=True)
#print(d)
#print(d[0])
d_out = 10*d.mean(axis=0)
np.savetxt(os.path.join('distance_' + out_dir, 'mean_distance_' + out_dir + UserInput.idx + '.dat'),d_out)
#print(d_out)
d_out.sort()
d_out=d_out[::-1]
#print(d_out)
