import numpy as np
import mdtraj as md
import sys, time, os
import argparse
import pandas as pd
from scipy.spatial.distance import euclidean
import glob


# Currently only working in python 2 due to MDAnalysis package.
# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description = 'Pull pdbs of clusters from time series', add_help=False) 

# List all possible user input
inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-top', action='store', dest='structure',help='Structure file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-traj', action='store', dest='trajectory',help='Trajectory',type=str,required=True)
inputs.add_argument('-sel', action='store', dest='sel', help='Atom selection',type=str,default='not element H')
inputs.add_argument('-i', action='store', dest='input',help='Text file containing cluster timeseries',type=str,required=True)
inputs.add_argument('-o', action='store', dest='outdir',help='File name for storing outdir',type=str,required=True)

UserInput=parser.parse_args()

# Map frames and clusters in dataframe
time_series = np.genfromtxt(UserInput.input)
num_frames = len(time_series)
labeled_traj = pd.DataFrame(columns = ['frame', 'cluster'])
labeled_traj['frame'] = np.arange(num_frames)
labeled_traj['cluster'] = time_series

# Load trajectory as matrix
topology = UserInput.structure
trajectory = UserInput.trajectory
t = md.load(trajectory,top=topology)
sel = t.topology.select(UserInput.sel)
t = t.atom_slice(sel)

# Format trajectory 
temp = t.xyz
frames = t.xyz.shape[0]
atoms = t.xyz.shape[1]
data = temp.reshape((frames,atoms*3))
data = data.astype('float64')
temp = []
t = []
centers = np.empty(int(max(time_series)) + 1)
# Restore all atoms for saving
t = md.load(trajectory,top=topology)


for i in range(0, int(max(time_series)) + 1):
    # We'll save each cluster's statistics in its own folder
    #directory = UserInput.outdir + '/cluster' + str(i)
    #os.makedirs(directory)
    cluster_string = labeled_traj.loc[labeled_traj['cluster'] == i].frame.values
    cluster_coords = data[cluster_string]
    mean = cluster_coords.mean(axis=0)
    distance = [euclidean(row,mean) for row in cluster_coords]
    rep = cluster_string[np.argmin(distance)]
    #t[rep].save_pdb(directory + '/rep.pdb')
    #t[cluster_string].save_pdb(directory + '/all.pdb')
    centers[i] = int(rep)
np.savetxt(UserInput.outdir + '/calculated_centers.txt', centers, fmt='%i')
