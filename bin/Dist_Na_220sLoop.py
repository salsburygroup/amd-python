#!/usr/bin/env python
import MDAnalysis
import argparse
import numpy as np
from scipy.spatial import distance
import hdbscan
import matplotlib.pyplot as plt
from Analysis.Cluster import Clusterer, Saver
from Analysis import TrajectoryReader
import os

print('Make sure to use the absolute path for output directory(-o)')
# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Compute closest mean distance between Na+ binding loop and Na+', add_help=False)

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

inputs.add_argument('-title',
                    action='store',
                    dest='title',
                    help='Title of the plot',
                    type=str,
                    default='Closest mean distance between Na$\mathregular{^+}$ and 220s loop')

inputs.add_argument('-tm',
                    action='store',
                    dest='timestep',
                    help='Timestep between two frames',
                    type=str,
                    default='0.1ns')

inputs.add_argument('-m',
                    action='store',
                    dest='method',
                    help='Clustering Method',
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

# Make output directory
if not os.path.exists(UserInput.out_dir):
    os.mkdir(UserInput.out_dir)

# Subtract the coordinate between the 
u = MDAnalysis.Universe(UserInput.structure, UserInput.trajectory)
sodiumLoop=u.select_atoms(UserInput.sel)
sodiumGrp=u.select_atoms('resname SOD')
T=len(u.trajectory)
data=np.zeros((T, len(sodiumLoop)*3))
data_no_noise=np.zeros((T, len(sodiumLoop)*3))
for i in range(T):
    u.trajectory[i] #Go to current frame
    Na_dist_mean=np.mean(distance.cdist(sodiumLoop.positions,sodiumGrp.positions,'euclidean'),axis=0)
    Na_min_index=np.argmin(Na_dist_mean)
    a=sodiumGrp.positions[Na_min_index]
    b=sodiumLoop.positions
    data[i]=(a-b).flatten()
    if np.min(Na_dist_mean)<12:
        data_no_noise[i]=(a-b).flatten()
np.savetxt(os.path.join(UserInput.out_dir, 'data.dat'), data, delimiter=' ')
np.savetxt(os.path.join(UserInput.out_dir, 'data_no_noise.dat'), data_no_noise, delimiter=' ')

# Cluster
#clusterer = hdbscan.HDBSCAN(min_cluster_size=10)
#clusterer = getattr(Clusterer, 'HDBSCAN')(data)
clusterer = getattr(Clusterer, UserInput.method)(data)
clusterer.fit()
#clusterer.labels
# Save Timeseries
Saver.TimeSeries(out_name=os.path.join(UserInput.out_dir, 'timeseries.txt'), labels=clusterer.labels).save()
# Plot time series
frames = np.arange(clusterer.labels.shape[0])
fig1=plt.figure(1,figsize=(8,4))
plt.scatter(frames, clusterer.labels, marker='+')
plt.xlabel('Frame(' + UserInput.timestep + ')')
plt.ylabel('Cluster')
plt.title(UserInput.title)
fig1.savefig(os.path.join(UserInput.out_dir, 'plot_timeseries.png'), dpi=200)
fig1.savefig(os.path.join(UserInput.out_dir, 'plot_timeseries.tiff'), dpi=600)
fig1.savefig(os.path.join(UserInput.out_dir, 'plot_timeseries.pdf'), bbox_inches='tight')

# Extract PDB and DCD
full_trajectory = TrajectoryReader.DCD(topology_path=UserInput.structure, trajectory_path=UserInput.trajectory).load()
        
Saver.PDB(out_name=os.path.join(UserInput.out_dir, 'clusters'), labels=clusterer.labels, trajectory=full_trajectory, atom_selection='all').save()
