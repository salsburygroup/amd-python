#! /usr/bin/env python

import glob
import numpy
import pandas
import re
import mdtraj
from scipy.spatial.distance import euclidean
from Analysis import AtomSelection, Correlation, Distance, Featurizer, Plotter, Saver, TrajectoryReader, TrajectoryProcessor

atom_selection = "name CA"

folders = glob.glob('*/stride*/')
for folder in folders:
    prefix = re.findall(r"(.*)/stride\d+", folder)[0]
    stride = re.findall(r"stride(\d+)", folder)[0]
    topology = glob.glob("{0}/*.pdb".format(prefix))[0]
    reps_dcd = glob.glob("{0}test{1}_reps.dcd".format(folder, stride))[0]
    stride_dcd = glob.glob("{0}test{1}.dcd".format(folder, stride))[0]
    kept_details_csv = glob.glob("{0}test{1}_kept_details.csv".format(folder, stride))[0]
    kept_time_series = glob.glob("{0}test{1}_kept_only_clustering.txt".format(folder, stride))[0]
    reps_trajectory = TrajectoryReader.DCD(
        topology_path=topology,
        trajectory_path=reps_dcd
    ).load()
    reps_trajectory = AtomSelection.Slice(trajectory=reps_trajectory, atom_selection=atom_selection).select()
    reps_trajectory = TrajectoryProcessor.Aligner(trajectory=reps_trajectory, atom_selection=atom_selection).process()
    kept_details = pandas.read_csv(kept_details_csv)
    kept_time_series = numpy.genfromtxt(kept_time_series)

    # Make minimum trajectory from kept clustering
    weights = kept_details['population'].astype(int)
    reps_original_trajectory = TrajectoryReader.DCD(
        topology_path=topology,
        trajectory_path=reps_dcd
    ).load()
    reps_original_trajectory = TrajectoryProcessor.Aligner(
            trajectory=reps_original_trajectory, atom_selection=atom_selection).process()
    num_frames = len(kept_time_series)
    labeled_traj = pandas.DataFrame(columns=['frame', 'cluster'])
    labeled_traj['frame'] = numpy.arange(num_frames)
    labeled_traj['cluster'] = kept_time_series
    final_weights = numpy.empty(int(max(kept_time_series)) + 1)
    data = Featurizer.XYZ(reps_trajectory).extract()

    with mdtraj.formats.DCDTrajectoryFile('{0}{1}_minimum_{2}_trajectory.dcd'.format(folder, prefix, stride),
                                          'w') as minimum_trajectory:
        for i in range(0, int(max(kept_time_series)) + 1):
            cluster_string = labeled_traj.loc[labeled_traj['cluster'] == i].frame.values
            cluster_weights = weights[cluster_string]
            final_weights[i] = numpy.sum(cluster_weights)
            cluster_coords = data[cluster_string]
            mean = cluster_coords.mean(axis=0)
            distance = [euclidean(row, mean) for row in cluster_coords]
            rep = cluster_string[numpy.argmin(distance)]
            minimum_trajectory.write(reps_original_trajectory.xyz[rep])
    numpy.savetxt('{0}{1}_minimum_{2}_weights.txt'.format(folder, prefix, stride), final_weights, fmt='%i')
    del reps_original_trajectory
    del reps_trajectory
    del data

    minimum_trajectory = TrajectoryReader.DCD(
        trajectory_path='{0}{1}_minimum_{2}_trajectory.dcd'.format(folder, prefix, stride),
        topology_path=topology
    ).load()
    minimum_trajectory = AtomSelection.Slice(trajectory=minimum_trajectory, atom_selection=atom_selection).select()
    minimum_trajectory = TrajectoryProcessor.Aligner(
        trajectory=minimum_trajectory, atom_selection=atom_selection).process()

    # Weighted correlation for minimu trajectory
    average = numpy.average(minimum_trajectory.xyz, axis=0, weights=final_weights)
    fluctuations = minimum_trajectory.xyz - average[numpy.newaxis, :]
    del average
    dots = numpy.zeros((minimum_trajectory.n_atoms, minimum_trajectory.n_atoms))
    for i in range(minimum_trajectory.n_frames):
        dot = numpy.dot(fluctuations[i, :, :], numpy.transpose(fluctuations[i, :, :]))
        dots = dots + dot * final_weights[i]
    del fluctuations
    dots = numpy.divide(dots, numpy.sum(final_weights))
    diagonal = numpy.diag(dots)
    normalization_matrix = numpy.outer(diagonal, diagonal)
    normalization_matrix = numpy.sqrt(normalization_matrix)
    reps_correlation_matrix = numpy.divide(dots, normalization_matrix)
    Plotter.UnityPColor(y=reps_correlation_matrix,
                        out_name="{0}{1}kept_{2}_correlation.png".format(folder, prefix, stride),
                        x_label=atom_selection,
                        y_label=atom_selection,
                        title="Kept Correlation Stride {0}".format(stride)
                        ).plot()
    Saver.Array(
        array=reps_correlation_matrix,
        out_name="{0}{1}_kept_{2}_correlation.txt".format(folder, prefix, stride)
    ).save()

    # Get strided trajectory
    stride_trajectory = TrajectoryReader.DCD(
        topology_path=topology,
        trajectory_path=stride_dcd
    ).load()
    stride_trajectory = AtomSelection.Slice(trajectory=stride_trajectory, atom_selection=atom_selection).select()
    stride_trajectory = TrajectoryProcessor.Aligner(
        trajectory=stride_trajectory, atom_selection=atom_selection).process()

    # Regular correlation for stride
    stride_correlation_matrix = Correlation.Pearson(trajectory=stride_trajectory).calculate()
    Plotter.UnityPColor(y=stride_correlation_matrix,
                        out_name="{0}{1}strided_{2}_correlation.png".format(folder, prefix, stride),
                        x_label=atom_selection,
                        y_label=atom_selection,
                        title="{0} Strided Correlation Stride {1}".format(prefix, stride)
                        ).plot()
    Saver.Array(
        array=stride_correlation_matrix,
        out_name="{0}{1}strided_{2}_correlation.txt".format(folder, prefix, stride)
    ).save()

    correlation_correlation = numpy.corrcoef(reps_correlation_matrix.ravel(), stride_correlation_matrix.ravel())[1][0]
    print("For {0}, stride {1}, correlation of correlation matrices is {2}".format(
        prefix, stride, correlation_correlation))

    # Weighted RMSF for kept
    average = numpy.average(minimum_trajectory.xyz, axis=0, weights=final_weights)
    fluctuations = minimum_trajectory.xyz - average[numpy.newaxis, :]
    sum_squares = numpy.sum(numpy.sum(numpy.square(fluctuations), axis=2)*final_weights[:, numpy.newaxis], axis=0)
    reps_rmsf = (sum_squares / numpy.sum(final_weights)) ** 0.5
    Plotter.Y(y=reps_rmsf,
              out_name="{0}{1}kept_{2}_rmsf.png".format(folder, prefix, stride),
              x_label=atom_selection,
              y_label='RMSF (nm)',
              title="{0} Kept RMSF Stride {1}".format(prefix, stride)
              ).plot()
    Saver.Array(array=reps_rmsf,
                out_name="{0}{1}kept_{2}_rmsf.txt".format(folder, prefix, stride)
                ).save()

    # Regular RMSF for strided
    stride_rmsf = Distance.RMSF(trajectory=stride_trajectory,
                                atom_selection=atom_selection
                                ).calculate()
    Plotter.Y(y=stride_rmsf,
              out_name="{0}{1}Stride_{2}_rmsf.png".format(folder, prefix, stride),
              x_label=atom_selection,
              y_label='RMSF (nm)',
              title="{0} Strided RMSF Stride {1}".format(prefix, stride)
              ).plot()

    Saver.Array(array=stride_rmsf,
                out_name="{0}{1}strided_{2}_rmsf.txt".format(folder, prefix, stride)
                ).save()
