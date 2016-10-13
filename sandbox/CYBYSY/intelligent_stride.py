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
import seaborn
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

# tracker will be our master record keeper, making sure we know where everything comes from
tracker = pandas.DataFrame(columns=['step', 'frame_in_rep_dcd', 'population'])
rep_trajectory_file = UserInput.output_prefix + '_reps.dcd'
step = 0
frame_counter = 0

with mdtraj.formats.DCDTrajectoryFile(
                UserInput.output_prefix + '.dcd', 'w') as strided_trajectory, mdtraj.formats.DCDTrajectoryFile(
                rep_trajectory_file, 'w') as rep_trajectory:
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
        step_dir = UserInput.output_prefix + '_step' + str(step)
        for i in range(0, int(max(cluster_labels)) + 1):
            cluster_string = labeled_traj.loc[labeled_traj['cluster'] == i].frame.values
            cluster_coords = data[cluster_string]
            mean = cluster_coords.mean(axis=0)
            distance = [euclidean(row, mean) for row in cluster_coords]
            rep = cluster_string[numpy.argmin(distance)]
            rep_trajectory.write(window[rep].xyz)
            tracker.loc[len(tracker)] = [step, frame_counter, len(cluster_coords)]
            frame_counter += 1
        del window_slice
        strided_trajectory.write(window[-1].xyz)
        del window
        step += 1

# Second round of clustering, keep track of populations along the way
reps_traj = mdtraj.load(rep_trajectory_file, top=UserInput.structure_file)
frames = reps_traj.n_frames
atoms = reps_traj.n_atoms
reps = reps_traj.xyz
reps = reps.reshape((frames, atoms * 3))
reps = reps.astype('float64')
clusterer = hdbscan.HDBSCAN(2)
tracker['final_label'] = clusterer.fit_predict(reps)
cluster_numbers = numpy.arange(start=-1, stop=max(tracker['final_label'] + 1))
population = numpy.empty(shape=(len(cluster_numbers), 1))
for index, cluster in enumerate(cluster_numbers):
    population[index] = int(sum(tracker[tracker.final_label == cluster].population))
norm_pop = population/sum(population)
relative_populations = numpy.column_stack((cluster_numbers, norm_pop))
numpy.save(UserInput.output_prefix + '_populations.txt', relative_populations)
for i in range(0, int(max(tracker['final_label'])) + 1):
    cluster_string = tracker.loc[tracker.final_label == i].frame_in_rep_dcd.values
    cluster_coords = reps[cluster_string.astype('int')]
    mean = cluster_coords.mean(axis=0)
    distance = [euclidean(row, mean) for row in cluster_coords]
    rep = cluster_string[numpy.argmin(distance)]
    reps_traj[rep].save_pdb(os.path.join(UserInput.output_prefix + '_rep' + str(i) + '.pdb'))

# Bar chart of populations
# width = 0.35
plt.bar(relative_populations[:,0], relative_populations[:,1])
plt.title('Relative populations')
plt.ylabel('Normalized population over original trajectory')
plt.xlabel('Cluster number')
plt.savefig(UserInput.output_prefix + '_rep_populations.png')

numpy.savetxt(UserInput.output_prefix + '_hdbscan_labels.txt', tracker.final_label.values)
tracker.to_csv(UserInput.output_prefix + '_details.csv')
# plt.xticks(relative_populations[:,0] + width/2., relative_populations[:,0])


