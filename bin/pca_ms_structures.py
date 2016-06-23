import numpy as np
import mdtraj
import sklearn.cluster
import sklearn.decomposition
import argparse
import matplotlib.pyplot as plt
from itertools import cycle
import sys
from scipy.spatial.distance import euclidean
from Analysis import AtomSelection, DimensionReduction, Featurizer, TrajectoryReader

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Find Gaussian wells in PCA space', add_help=False)

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
inputs.add_argument('-sel',
                    action='store',
                    dest='sel',
                    help='Atom selection',
                    type=str,
                    default='name CA'
                    )
inputs.add_argument('-n',
                    action='store',
                    dest='max_components',
                    help='Number of components to project and plot in 2D',
                    type=int,
                    default=2
                    )
inputs.add_argument('-o',
                    action='store',
                    dest='out_name',
                    help='Output folder',
                    type=str,
                    required=True
                    )

# Parse into useful form
UserInput = parser.parse_args()

# Process trajectory
trajectory = TrajectoryReader.DCD(topology_path=UserInput.structure, trajectory_path=UserInput.trajectory).load()
trajectory = AtomSelection.Slice(trajectory=trajectory, atom_selection=UserInput.sel).select()
assert isinstance(trajectory, mdtraj.Trajectory)
xyz = Featurizer.XYZ(trajectory=trajectory).extract()

projection, components, explained_variance = DimensionReduction.PCA(coordinates=xyz).reduce()
for i in range(0, UserInput.max_components-1):
    for j in range(i+1, UserInput.max_components):
        clusterer = sklearn.cluster.MeanShift(n_jobs=-1, cluster_all=True)
        temp = np.column_stack((projection[:, i], projection[:, j]))
        labels = clusterer.fit_predict(temp)
        cluster_centers = clusterer.cluster_centers_
        labels_unique = np.unique(labels)
        n_clusters_ = max(labels_unique)+1
        center_indices = []
        for center in cluster_centers:
            distance = [euclidean(row, center) for row in temp]
            center_indices.append(labels[np.argmin(distance)])
        trajectory[center_indices].save(UserInput.out_name + '/MScenter_PCA_{0}_{1}.pdb'.format(i, j))
        colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
        for k, col in zip(range(n_clusters_), colors):
            my_members = labels == k
            cluster_center = cluster_centers[k]
            plt.plot(projection[my_members, i], projection[my_members, j], col + '.')
            plt.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
                     markeredgecolor='k', markersize=14)
            
        progress = "\r Calculations running for PCA {0} of {1}".format(i + 1, UserInput.max_components)  # status
        sys.stdout.write(progress)
        sys.stdout.flush()  # report status to terminal output
        plt.title('Estimated number of clusters: %d' % n_clusters_)
        plt.savefig(UserInput.out_name + '/PCAMeanShift_{0}_{1}.png'.format(i, j))
        plt.clf()

print("\rColor order for matching pdbs to centers is bgrcmykbgrcmykbgrcmykbgrcmyk")
