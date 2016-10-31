#!/usr/bin/env python

# Ideally, you'll use the same atom selection you would for a correlation matrix.
# I suggest 'name CA' for proteins and 'not element H' for nucleic acids.

from Analysis import AtomSelection, Cluster, Correlation, Saver, TrajectoryReader
import argparse
import matplotlib.pyplot as plt
# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Calculate Propagating Correlation Matrix', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-top', action='store', dest='structure',
                    help='Structure file corresponding to trajectory', type=str, required=True)
inputs.add_argument('-traj', action='store', dest='trajectory', help='Trajectory', type=str, required=True)
inputs.add_argument('-sel', action='store', dest='sel', help='Atom selection', type=str, default='not element H')
inputs.add_argument('-tau', action='store', dest='tau', help='lag time', type=int, default=1)
inputs.add_argument('-o', action='store', dest='out_name', help='Output prefix', type=str, required=True)
inputs.add_argument('-c',
                    action='store_true',
                    dest='also_cluster',
                    help='Perform correlation clustering as well.',
                    )
inputs.add_argument('-d',
                    action='store_true',
                    dest='correlation_as_distance',
                    help='Use 1 minus abs(correlation) matrix as distance for clustering.',
                    )
inputs.add_argument('-p',
                    action='store_true',
                    dest='dot_products_as_distance',
                    help='Use dot products as similarity matrix for clustering.',
                    )
inputs.add_argument('-m',
                    action='store',
                    dest='min_cluster_size',
                    help='Minimum cluster and neighborhood size for HDBSCAN clustering.',
                    type=int
                    )

# Parse into useful form
UserInput = parser.parse_args()
trajectory = TrajectoryReader.DCD(topology_path=UserInput.structure, trajectory_path=UserInput.trajectory).load()
trajectory = AtomSelection.Slice(trajectory=trajectory, atom_selection=UserInput.sel).select()
# Execute calculation
cp = Correlation.Propagator(trajectory=trajectory, tau=UserInput.tau)
average_dot, average_delta, dot_average_delta = cp.calculate()

# Save text results
Saver.Array(array=average_dot, out_name=UserInput.out_name + '_average_dot.txt').save()
Saver.Array(array=average_delta, out_name=UserInput.out_name + '_average_delta.txt').save()
Saver.Array(array=dot_average_delta, out_name=UserInput.out_name + '_dot_average_delta.txt').save()

# Save pretty pictures
# need to add matshow to Plotter
plt.matshow(average_dot)
plt.xlabel('Atom')
plt.ylabel('Atom')
plt.title('Average dot product')
plt.colorbar()
plt.savefig(UserInput.out_name + '_average_dot.png')
plt.close()

#plt.figure()
#plt.matshow(dot_average_delta)
#plt.xlabel('Atom')
#plt.ylabel('Atom')
#plt.title('Dot product of average delta')
#plt.colorbar(format='%.0e')
#plt.savefig(UserInput.out_name + '_dot_average_delta.png')
#plt.close()

if UserInput.also_cluster:
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

    # temporarily having this do all three for deciding which makes the most sense
    labels = Correlation.Clustering.cluster(correlation_matrix=average_dot,  input_type='correlation')
    Cluster.Saver.ClusterFrames(out_name=UserInput.out_name + '_residue_groups_correlation.txt', labels=labels).save()
    Correlation.Clustering.visualize(
        labels=labels, pdb_file=UserInput.structure, out_name=UserInput.out_name + '_correlation_render.tga'
    )
    labels = Correlation.Clustering.cluster(correlation_matrix=average_dot, input_type='similarity')
    Cluster.Saver.ClusterFrames(out_name=UserInput.out_name + '_residue_groups_similarity.txt', labels=labels).save()
    Correlation.Clustering.visualize(
        labels=labels, pdb_file=UserInput.structure, out_name=UserInput.out_name + '_similarity_render.tga'
    )
    labels = Correlation.Clustering.cluster(correlation_matrix=average_dot, input_type='dots')
    Cluster.Saver.ClusterFrames(out_name=UserInput.out_name + '_residue_groups_similarity.txt', labels=labels).save()
    Correlation.Clustering.visualize(
        labels=labels, pdb_file=UserInput.structure, out_name=UserInput.out_name + '_dots_render.tga'
    )