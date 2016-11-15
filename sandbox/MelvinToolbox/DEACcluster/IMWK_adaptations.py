#! /usr/bin/env python

# Script for Rescaled IMWK-Means followed by X, where X is a clustering method.


import argparse
import os
from Analysis import TrajectoryReader, Featurizer, AtomSelection
from Analysis.Cluster import Clusterer, Plotter, Saver, Scorer
from sklearn.cluster import MiniBatchKMeans
from sklearn import mixture
from sklearn import cluster


# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
    description='Runs cluster trials over variety of methods and metrics', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument(
    '-top',
    action='store',
    dest='structure',
    help='Structure file corresponding to trajectory',
    type=str,
    required=True)
inputs.add_argument('-traj', action='store', dest='trajectory', help='Trajectory', type=str, required=True)
inputs.add_argument('-sel', action='store', dest='sel', help='atoms', type=str, required=True)
inputs.add_argument('-o', action='store', dest='out_name', help='Output directory', type=str, required=True)

# Parse into useful form
UserInput = parser.parse_args()

trajectory = TrajectoryReader.DCD(trajectory_path=UserInput.trajectory, topology_path=UserInput.structure).load()
trajectory = AtomSelection.Slice(trajectory=trajectory, atom_selection=atom_selection).select()
trajectory_2d = Featurizer.XYZ(trajectory=trajectory).extract()


def make_directory(method_name):
    output_path = os.path.join(UserInput.out_name, method_name)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    return output_path


def finishing(trial_labels, cwd):
    Saver.TimeSeries(out_name=os.path.join(cwd, 'timeseries.txt'), labels=labels).save()
    Plotter.TimeSeries(out_name=os.path.join(cwd, 'timeseries.png'), labels=labels).plot()
    Saver.ClusterFrames(out_name=os.path.join(cwd, 'row_format.txt'), labels=trial_labels).save()
    silhouette_score = Scorer.Silhouette(labels=trial_labels, data=trajectory_2d).evaluate()
    Saver.Score(out_name=os.path.join(cwd, 'silhouette.txt'), score=silhouette_score).save()
    Saver.PDB(
        out_name=os.path.join(cwd, 'clusters'),
        labels=trial_labels, trajectory=trajectory, atom_selection='all'
    ).save()

# Rescaled IMWK-Means data
IMWK_labels, IMWK_center, IMWK_data, IMWK_weights, IMWK_optimal_k = Clusterer.IMWKRescaled(
    trajectory_2d=trajectory_2d).fit()

# K-Means (A-H)
directory = make_directory('AmorimHennig')
try:
    labels = MiniBatchKMeans(n_clusters=IMWK_optimal_k, n_init=5).fit_predict(IMWK_data)
    finishing(labels, directory)
except Exception:
    pass

# HDBSCAN
directory = make_directory('HDBSCAN')
try:
    labels = Clusterer.HDBSCAN(trajectory_2d=IMWK_data).fit()
    finishing(labels, directory)
except Exception:
    pass

# # GMM
# directory = make_directory('GMM')
# labels = mixture.GMM(n_components=IMWK_optimal_k, covariance_type='tied').fit_predict(IMWK_data)
# finishing(labels, directory)
#
# # VBGMM
# directory = make_directory('VBGMM')
# labels = mixture.VBGMM(n_components=IMWK_optimal_k, covariance_type='tied').fit_predict(IMWK_data)
# finishing(labels, directory)

# MeanShift
directory = make_directory('MeanShift')
try:
    labels = cluster.MeanShift(n_jobs=-1).fit_predict(IMWK_data)
    finishing(labels, directory)
except Exception:
    pass

# AffinityPropagation
directory = make_directory('AffinityPropagation')
try:
    labels = cluster.AffinityPropagation(damping=0.5).fit_predict(IMWK_data)
    finishing(labels, directory)
except Exception:
    pass
