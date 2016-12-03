#!/usr/env/ Python

# Ryan Melvin


import hdbscan
import os
import pandas
import numpy
import mdtraj
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
from scipy.spatial.distance import euclidean
import time

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

start_time = time.time()

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
strided_trajectory = mdtraj.load(UserInput.output_prefix + '.dcd', top=UserInput.structure_file)
reps_trajectory = mdtraj.load(rep_trajectory_file, top=UserInput.structure_file)
strided_trajectory = strided_trajectory.superpose(strided_trajectory)
reps_trajectory = reps_trajectory.superpose(strided_trajectory)
strided_frames = strided_trajectory.n_frames
atoms = strided_trajectory.n_atoms
strided_trajectory = strided_trajectory.xyz
strided_trajectory = strided_trajectory.reshape((strided_frames, atoms * 3))
rep_frames = reps_trajectory.n_frames
reps = reps_trajectory.xyz
reps = reps.reshape((rep_frames, atoms * 3))
reps = reps.astype('float64')
both_trajectories = numpy.row_stack((strided_trajectory, reps))
clusterer = hdbscan.HDBSCAN(min_cluster_size=2)
labels = clusterer.fit_predict(both_trajectories)
strided_labels = labels[0:strided_frames-1]
rep_labels = labels[strided_frames-1:-1]
rep_labeled_traj = pandas.DataFrame(columns=['frame', 'cluster'])
rep_labeled_traj['frame'] = numpy.arange(len(rep_labels))
rep_labeled_traj['cluster'] = rep_labels
end_time = time.time()
print("Elapsed time was %g seconds" % (end_time - start_time))
strided_unique_clusters = [cluster_number for cluster_number in strided_labels if cluster_number not in rep_labels]
strided_unique_clusters = numpy.unique(strided_unique_clusters)
strided_noise = [ind for ind, cluster in enumerate(strided_labels) if cluster == -1]
rep_unique_clusters = [cluster_number for cluster_number in rep_labels if cluster_number not in strided_labels]
rep_unique_clusters = numpy.unique(rep_unique_clusters)
rep_noise = [ind for ind, cluster in enumerate(rep_labels) if cluster == -1]
print(str(len(rep_unique_clusters) + len(rep_noise)) + ' things missed in striding\n'
      + '\t' + str(len(rep_unique_clusters)) + ' are from labeled clusters\n'
      + '\t' + str(len(rep_noise)) + ' are noise\n')
print(str(len(strided_unique_clusters) + len(strided_noise)) + ' things only in strided trajectory\n'
      + '\t' + str(len(strided_unique_clusters)) + ' are from labeled clusters\n'
      + '\t' + str(len(strided_noise)) + ' are noise\n')

matplotlib.pyplot.figure()
matplotlib.pyplot.scatter(numpy.arange(len(labels)), labels, marker='+')
matplotlib.pyplot.xlabel('Frame (only meaningful to line)')
matplotlib.pyplot.ylabel('Cluster')
matplotlib.pyplot.title('Intelligent Stride Final Round')
matplotlib.pyplot.savefig(UserInput.output_prefix + '_final_round_clusters.png')
matplotlib.pyplot.close()
tracker['final_label'] = rep_labels
if not os.path.exists(UserInput.output_prefix + '_MissedStructures'):
    os.mkdir(UserInput.output_prefix + '_MissedStructures')
for label in rep_unique_clusters:
    cluster_string = rep_labeled_traj.loc[rep_labeled_traj['cluster'] == label].frame.values
    cluster_coords = reps[cluster_string]
    mean = cluster_coords.mean(axis=0)
    distance = [euclidean(row, mean) for row in cluster_coords]
    rep = cluster_string[numpy.argmin(distance)]
    population = int(sum(tracker[tracker.final_label == label].population))
    reps_trajectory[rep].save_pdb(UserInput.output_prefix + '_MissedStructures/' + '/rep_' + str(label)
                                  + '_population_' + str(population) + '.pdb')
for frame_number in rep_noise:
    population = int(sum(tracker[tracker.frame_in_rep_dcd == frame_number].population))
    reps_trajectory[frame_number].save_pdb(UserInput.output_prefix + '_MissedStructures/' + '/noise_frame'
                                           + str(frame_number) + '_population_' + str(population) + '.pdb')

tracker.to_csv(UserInput.output_prefix + '_kept_details.csv')
strided_clusters, strided_cluster_counts = numpy.unique(strided_labels, return_counts=True)
stride_details = pandas.DataFrame(columns=['cluster', 'population'])
stride_details['cluster'] = strided_clusters
stride_details['population'] = strided_cluster_counts
stride_details.to_csv(UserInput.output_prefix + '_stride_details.csv')

summary = pandas.DataFrame(columns=['clusters_only_in_kept',
                                    'fraction_of_original',
                                    'clusters_only_in_strided',
                                    'fraction_of_strided',
                                    'clusters_in_both'])
clusters_only_in_kept = len(rep_unique_clusters)
total_frames = strided_frames * UserInput.stride
fraction_original  = sum(tracker[tracker['final_label'].isin(rep_unique_clusters)]['population'])/total_frames
clusters_only_in_strided = len(strided_unique_clusters)
fraction_strided = \
    sum(stride_details[stride_details['cluster'].isin(strided_unique_clusters)]['population'])/strided_frames
clusters_in_both = [cluster_number for cluster_number in strided_labels if cluster_number in rep_labels]
clusters_in_both = len(clusters_in_both)
summary.loc[0] = [clusters_only_in_kept, fraction_original, clusters_only_in_strided,
                                    fraction_strided, clusters_in_both]
summary.to_csv(UserInput.output_prefix + '_summary.csv')