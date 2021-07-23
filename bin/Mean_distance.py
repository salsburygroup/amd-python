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
print(d[0])
d_out = 10*d.mean(axis=0)
np.savetxt(os.path.join('distance_' + out_dir, 'mean_distance_' + out_dir + '.dat'),d_out)
#print(d_out)
d_out.sort()
d_out=d_out[::-1]
#print(d_out)

# Plot
fig1=plt.figure(1,figsize=(16,4))
plt.plot(d_out)
plt.xlabel('residue pair number')
plt.ylabel('Distance (Å)')
plt.title(UserInput.title)
fig1.savefig(os.path.join('distance_' + out_dir, 'mean_distance_' + out_dir + '.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
fig1.savefig(os.path.join('distance_' + out_dir, 'mean_distance_' + out_dir + '.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
fig1.savefig(os.path.join('distance_' + out_dir, 'mean_distance_' + out_dir + '.pdf'), bbox_inches='tight')

# Histogram
num_of_bins=int(1+np.log(len(d_out))/np.log(2))
#print(num_of_bins)
counts, bins=np.histogram(d_out,num_of_bins)
#print(counts)

# Plot histogram
fig2=plt.figure(2,figsize=(16,4))
plt.hist(bins[:-1],bins,weights=counts/len(d_out))
plt.xlabel('Distance (Å)')
plt.ylabel('Frequency')
plt.title(UserInput.title)
fig2.savefig(os.path.join('distance_' + out_dir, 'hist_mean_distance_' + out_dir + '.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
fig2.savefig(os.path.join('distance_' + out_dir, 'hist_mean_distance_' + out_dir + '.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
fig2.savefig(os.path.join('distance_' + out_dir, 'hist_mean_distance_' + out_dir + '.pdf'), bbox_inches='tight')
