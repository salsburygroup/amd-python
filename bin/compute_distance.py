#! /usr/bin/env python

import argparse
import os
import numpy
import mdtraj as md
import pandas as pd
from scipy.spatial import distance
#from Analysis import Plotter, Saver

# Jiajie Xiao
# 07.19.2017

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
    description='Compute (min,max,mean or all pairwise) distance of selected  pairs of atoms for given MD trajectories', add_help=False
)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-top',
                    action='store',
                    dest='structure',
                    help='Structure file corresponding to trajectory',
                    type=str,
                    required=True
                    )
inputs.add_argument('-traj',
                    action='store',
                    dest='trajectory',
                    help='Trajectory',
                    type=str,
                    required=True
                    )
inputs.add_argument('-stride',
                    action='store',
                    dest='stride',
                    help='stride to use',
                    type=int,
                    default=1
                    )
inputs.add_argument('-periodic',
                    action='store',
                    dest='periodicTrue',
                    help='periodic true of false',
                    type=bool,
                    default=True
                    )
inputs.add_argument('-sel1',
                    action='store',
                    dest='sel1',
                    help='selection 1',
                    type=str,
                    required=True
                    )
inputs.add_argument('-sel2',
                    action='store',
                    dest='sel2',
                    help='selection 2',
                    type=str,
                    required=True
                    )
inputs.add_argument('-mode',
                    action='store',
                    dest='mode',
                    help='output mode: min, mean, max, verbose, com, cog',
                    type=str,
                    default='verbose'
                    )
inputs.add_argument('-o',
                    action='store',
                    dest='out_name',
                    help='Output folder for data',
                    type=str,
                    required=True
                    )
# Parse into useful form
UserInput = parser.parse_args()
traj = md.load(UserInput.trajectory, top=UserInput.structure, stride=UserInput.stride)
atom_pairs = traj.topology.select_pairs(selection1=UserInput.sel1,selection2=UserInput.sel2)

if UserInput.mode == 'com':
    num_frames = traj.n_frames
    traj_sel1 = traj.atom_slice(traj.topology.select(UserInput.sel1))
    center_of_mass = md.compute_center_of_mass(traj_sel1)
    coor_sel2 = traj.atom_slice(traj.topology.select(UserInput.sel2)).xyz
    d_out = 10*numpy.array([numpy.min(distance.cdist(numpy.array([center_of_mass[frame]]), coor_sel2[frame],'euclidean')) for frame in range(num_frames)])
    numpy.savetxt('distance_COM'+UserInput.out_name+'_min.dat',d_out)
elif UserInput.mode == 'cog':
    num_frames = traj.n_frames
    traj_sel1 = traj.atom_slice(traj.topology.select(UserInput.sel1))
    center_geometry = traj_sel1.xyz.mean(axis=1)
    coor_sel2 = traj.atom_slice(traj.topology.select(UserInput.sel2)).xyz
    d_out = 10 * numpy.array([numpy.min(distance.cdist(numpy.array([center_geometry[frame]]), coor_sel2[frame], 'euclidean')) for frame in range(num_frames)])
    numpy.savetxt('distance_COG'+UserInput.out_name+'_min.dat',d_out)
else:
    d = md.compute_distances(traj,atom_pairs,periodic=UserInput.periodicTrue,opt=True)
    if UserInput.mode == 'min':
        d_out = 10*d.min(axis=1)
        numpy.savetxt('distance_'+UserInput.out_name+'_min.dat',d_out)
    elif UserInput.mode == 'max':
        d_out = 10*d.max(axis=1)
        numpy.savetxt('distance_'+UserInput.out_name+'_max.dat',d_out)
    elif UserInput.mode == 'mean':
        d_out = 10*d.mean(axis=1)
        numpy.savetxt('distance_'+UserInput.out_name+'_mean.dat',d_out)
    else:
        top, bonds = traj.top.to_dataframe()
        id0 = top.ix[atom_pairs[:,0]].resSeq
        id1 = top.ix[atom_pairs[:,1]].resSeq
        name0 = top.ix[atom_pairs[:,0]].name
        name1 = top.ix[atom_pairs[:,1]].name

        columns = pd.MultiIndex.from_arrays([id0, name0, id1, name1])
        distance = pd.DataFrame(10*d,columns=columns)
        distance.to_csv('distance_'+UserInput.out_name+'_verbose.dat')