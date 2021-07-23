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

inputs.add_argument('-m',
                    action='store',
                    dest='mode',
                    help='output mode: min, mean_atoms, mean_steps, max, verbose, com, cog',
                    type=str,
                    default='verbose')

inputs.add_argument('-tm',
                    action='store',
                    dest='time',
                    help='time interval(ps)',
                    type=float,
                    default=10)

inputs.add_argument('-title',
                    action='store',
                    dest='title',
                    help='Title of the plot',
                    type=str,
                    default='Mean distance between two dimers')

inputs.add_argument('-stride',
                    action='store',
                    dest='stride',
                    help='stride to use',
                    type=int,
                    required=True)

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output folder for data',
                    type=str,
                    required=True)

# Parse into useful form
UserInput = parser.parse_args()
traj = md.load(UserInput.trajectory, top=UserInput.structure, stride=UserInput.stride)
atom_pairs = traj.topology.select_pairs(selection1=UserInput.sel1,selection2=UserInput.sel2)

# Make output directory
out_dir=UserInput.out_dir
fname = './' + 'distance_' + out_dir
if not os.path.exists(fname):
    os.mkdir(fname)

# Make time axis
t=0.001*UserInput.time*UserInput.stride*np.array(range(1,traj.n_frames+1))

# All kinds of modes
if UserInput.mode == 'com':
    num_frames = traj.n_frames
    traj_sel1 = traj.atom_slice(traj.topology.select(UserInput.sel1))
    center_of_mass1 = md.compute_center_of_mass(traj_sel1)
    traj_sel2 = traj.atom_slice(traj.topology.select(UserInput.sel2))
    center_of_mass2 = md.compute_center_of_mass(traj_sel2)
    d_out = 10*np.array([distance.cdist(np.array([center_of_mass1[frame]]), np.array([center_of_mass2[frame]]), 'euclidean')[0] for frame in range(num_frames)])
    np.savetxt(os.path.join('distance_' + out_dir, 'COM_distance_' + out_dir + '.dat'),d_out)

elif UserInput.mode == 'com_min':
    num_frames = traj.n_frames
    traj_sel1 = traj.atom_slice(traj.topology.select(UserInput.sel1))
    center_of_mass = md.compute_center_of_mass(traj_sel1)
    coor_sel2 = traj.atom_slice(traj.topology.select(UserInput.sel2)).xyz
    d_out = 10*np.array([np.min(distance.cdist(np.array([center_of_mass[frame]]), coor_sel2[frame],'euclidean')) for frame in range(num_frames)])
    np.savetxt(os.path.join('distance_' + out_dir, 'COM_min_distance_' + out_dir + '.dat'),d_out)

elif UserInput.mode == 'cog':
    num_frames = traj.n_frames
    traj_sel1 = traj.atom_slice(traj.topology.select(UserInput.sel1))
    center_geometry = traj_sel1.xyz.mean(axis=1)
    coor_sel2 = traj.atom_slice(traj.topology.select(UserInput.sel2)).xyz
    d_out = 10*np.array([np.min(distance.cdist(np.array([center_geometry[frame]]), coor_sel2[frame], 'euclidean')) for frame in range(num_frames)])
    np.savetxt(os.path.join('distance_' + out_dir, 'COG_min_distance_' + out_dir + '.dat'),d_out)

elif UserInput.mode == 'verbose':
    top, bonds = traj.top.to_dataframe()
    id0 = top.ix[atom_pairs[:,0]].resSeq
    id1 = top.ix[atom_pairs[:,1]].resSeq
    name0 = top.ix[atom_pairs[:,0]].name
    name1 = top.ix[atom_pairs[:,1]].name
    columns = pd.MultiIndex.from_arrays([id0, name0, id1, name1])
    distance = pd.DataFrame(10*d,columns=columns)
    distance.to_csv('distance_'+out_dir+'_verbose.dat')

else:
    d = md.compute_distances(traj,atom_pairs,periodic=UserInput.periodic,opt=True)
    if UserInput.mode == 'min':
        d_out = 10*d.min(axis=1)
        np.savetxt(os.path.join('distance_' + out_dir, 'Min_distance_' + out_dir + '.dat'),d_out)

    elif UserInput.mode == 'max':
        d_out = 10*d.max(axis=1)
        np.savetxt(os.path.join('distance_' + out_dir, 'Max_distance_' + out_dir + '.dat'),d_out)

    elif UserInput.mode == 'mean_atoms':
        d_out = 10*d.mean(axis=1)
        np.savetxt(os.path.join('distance_' + out_dir, 'Mean_atoms_distance_' + out_dir + '.dat'),d_out)
    
    elif UserInput.mode == 'mean_steps':
        d_out = 10*d.mean(axis=0)
        np.savetxt(os.path.join('distance_' + out_dir, 'Mean_steps_distance_' + out_dir + '.dat'),d_out)
        d_out.sort()
        d_out=d_out[::-1]

#Xlabel
if UserInput.mode == 'mean_steps':
    Xlabel = 'residue pair number'
else:
    Xlabel = 'time (ns)'

# Plot
fig1=plt.figure(1,figsize=(16,4))
plt.plot(t, d_out, alpha=0.9, linewidth=0.5)
plt.xlabel(Xlabel)
plt.ylabel('Distance (Å)')
plt.title(UserInput.title)
fig1.savefig(os.path.join('distance_' + out_dir, UserInput.mode + '_distance_' + out_dir + '.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
fig1.savefig(os.path.join('distance_' + out_dir, UserInput.mode + '_distance_' + out_dir + '.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
fig1.savefig(os.path.join('distance_' + out_dir, UserInput.mode + '_distance_' + out_dir + '.pdf'), bbox_inches='tight')
