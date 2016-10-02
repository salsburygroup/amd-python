#!/usr/env/ Python

# Ryan Melvin

import os
import hdbscan
import pandas
import numpy
import mdtraj
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial.distance import euclidean


# Initialize parser for user input
parser = argparse.ArgumentParser(description='Run intelligently sampled ACEMD using HDBSCAN', add_help=False)
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-t',
                    '--trajectory',
                    action='store',
                    dest='trajectory_file',
                    help='Trajectory to be strided',
                    type=str,
                    required=True)
inputs.add_argument('-s',
                    '--structure',
                    action='store',
                    dest='structure_file',
                    help='Structure file corresponding to trajectory',
                    type=str,
                    required=True)
inputs.add_argument('-o',
                    '--output_prefix',
                    action='store',
                    dest='output_prefix',
                    help='Output prefix',
                    type=str,
                    required=True)
inputs.add_argument('-sel',
                    action='store',
                    dest='atom_selection',
                    help='Atom selection',
                    type=str,
                    default='not water')
inputs.add_argument('-stride',
                    action='store',
                    dest='stride',
                    help='Rate of resampling',
                    type=int,
                    default=100)
inputs.add_argument('-min',
                    action='store',
                    dest='minimum_cluster_size',
                    help='Minimum size of an HDBSCAN cluster. Default is 2, but higher is more efficient.',
                    type=int,
                    default=2)
UserInput = parser.parse_args()

step = 0
with mdtraj.formats.DCDTrajectoryFile(UserInput.output_prefix + '.dcd', 'w') as strided_trajectory:
    for window in mdtraj.iterload(UserInput.trajectory_file, top=UserInput.structure_file, chunk=UserInput.stride):
        window_slice = window.atom_slice(window.topology.select(UserInput.atom_selection))
        window_slice = window_slice.superpose(window_slice)
        temp = window_slice.xyz
        frames = window_slice.n_frames
        atoms = window_slice.n_atoms
        data = temp.reshape((frames, atoms * 3))
        data = data.astype('float64')
        del temp
        clusterer = hdbscan.HDBSCAN(UserInput.minimum_cluster_size)
        cluster_labels = clusterer.fit_predict(data)
        num_frames = len(cluster_labels)
        labeled_traj = pandas.DataFrame(columns=['frame', 'cluster'])
        labeled_traj['frame'] = numpy.arange(num_frames)
        labeled_traj['cluster'] = cluster_labels
        step_dir = UserInput.output_prefix + 'step' + str(step)
        if not os.path.exists(step_dir):
            os.mkdir(step_dir)
        for i in range(0, int(max(cluster_labels)) + 1):
            cluster_string = labeled_traj.loc[labeled_traj['cluster'] == i].frame.values
            cluster_coords = data[cluster_string]
            mean = cluster_coords.mean(axis=0)
            distance = [euclidean(row, mean) for row in cluster_coords]
            rep = cluster_string[numpy.argmin(distance)]
            window[rep].save_pdb(os.path.join(step_dir, 'rep' + str(i) + '.pdb'))
        plt.figure()
        plt.scatter(numpy.arange(frames), cluster_labels, marker='+')
        plt.xlabel('Frame')
        plt.ylabel('Cluster')
        plt.title('HDBSCAN')
        plt.savefig(os.path.join(step_dir, 'hdbscan_timeseries.png'))
        plt.close()
        del window_slice
        strided_trajectory.write(window[-1].xyz)
        del window
        step += 1
