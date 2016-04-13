import numpy as np
import mdtraj as md
import sklearn.cluster
import sklearn.decomposition
import argparse
import matplotlib.pyplot as plt
from itertools import cycle

# Currently only working in python 2 due to MDAnalysis package.
# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description = 'Run and score hdbscan clustering', add_help=False) 

# List all possible user input
inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-top', action='store', dest='structure',help='Structure file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-traj', action='store', dest='trajectory',help='Trajectory',type=str,required=True)
inputs.add_argument('-sel', action='store', dest='sel', help='Atom selection',type=str,default='name CA')
inputs.add_argument('-min', action='store', dest='min', help='minimum cluster membership',type=int,default=10)
inputs.add_argument('-o', action='store', dest='out_name',help='Output file',type=str,required=True)

# Parse into useful form
UserInput=parser.parse_args()

topology = UserInput.structure
trajectory = UserInput.trajectory
t = md.load(trajectory,top=topology)
sel = t.topology.select(UserInput.sel)
t = t.atom_slice(sel)

pca=sklearn.decomposition.PCA(n_components=2)
reduced_cartesian = pca.fit_transform(t.xyz.reshape(t.n_frames, t.n_atoms * 3))

bins = int(np.round(np.log2(t.xyz.shape[0])) + 1)
H, xedges, yedges = np.histogram2d(reduced_cartesian[:,0], reduced_cartesian[:,1], bins=bins)
H[H==0] = np.nan
E = -0.6 * np.log(H)
Em = np.ma.masked_where(np.isnan(E), E)

plt.pcolormesh(xedges, yedges, Em.T)
plt.savefig(UserInput.out_name + '/free_energy.png')
plt.close()

clusterer = sklearn.cluster.MeanShift()
labels = clusterer.fit_predict(reduced_cartesian)
cluster_centers = clusterer.cluster_centers_
labels_unique = np.unique(labels)
n_clusters_ = len(labels_unique)
print("number of estimated clusters : %d" % n_clusters_)

plt.figure()
plt.clf()

colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
for k, col in zip(range(n_clusters_), colors):
    my_members = labels == k
    cluster_center = cluster_centers[k]
    plt.plot(reduced_cartesian[my_members, 0], reduced_cartesian[my_members, 1], col + '.')
    plt.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=14)
plt.title('Estimated number of clusters: %d' % n_clusters_)
plt.savefig(UserInput.out_name + '/MeanShift.png')

