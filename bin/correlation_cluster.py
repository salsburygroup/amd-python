#!/usr/bin/env python

# Currently works with protein CA only

import argparse
from Analysis import AtomSelection, Correlation, TrajectoryReader
from Analysis.Cluster import Saver


# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Perform correlation clustering and generate images', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-pdb',
                    action='store',
                    dest='structure',
                    help='PDB corresponding to trajectory. MUST be PDB for this script.',
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
inputs.add_argument('-sel',
                    action='store',
                    dest='sel',
                    help='Atom selection. Must be name CA currently.',
                    type=str,
                    default='name CA'
                    )
inputs.add_argument('-m',
                    action='store',
                    dest='min_cluster_size',
                    help='Minimum cluster and neighborhood size for HDBSCAN',
                    type=int
                    )
inputs.add_argument('-d',
                    action='store_true',
                    dest='correlation_as_distance',
                    help='Use 1 minus abs(correlation) matrix as distance',
                    )
inputs.add_argument('-p',
                    action='store_true',
                    dest='dot_products_as_distance',
                    help='Use dot products as similarity matrix.',
                    )
inputs.add_argument('-o',
                    action='store',
                    dest='out_name',
                    help='Output prefix for vmd files and image',
                    type=str,
                    required=True
                    )

# Parse into useful form
UserInput = parser.parse_args()

# Process trajectory
trajectory = TrajectoryReader.DCD(topology_path=UserInput.structure, trajectory_path=UserInput.trajectory).load()
trajectory = AtomSelection.Slice(trajectory=trajectory, atom_selection=UserInput.sel).select()

if UserInput.correlation_as_distance:
    input_type = 'similarity'
elif UserInput.dot_products_as_distance:
    input_type = 'dots'
else:
    input_type = 'correlation'

if UserInput.min_cluster_size:
    min_membership = UserInput.min_cluster_size
else:
    min_membership = None

labels = Correlation.Clustering(
    trajectory=trajectory, input_type=input_type, minimum_membership=min_membership
    ).calculate()
Saver.ClusterFrames(out_name=UserInput.out_name + '_residue_groups.txt', labels=labels).save()
Correlation.Clustering.visualize(labels=labels, pdb_file=UserInput.structure, out_name=UserInput.out_name)
