import os
import numpy as np
import mdtraj as md
import hdbscan
from sklearn.metrics import silhouette_samples, silhouette_score
import matplotlib
matplotlib.use('Agg') #For use on DEAC cluster
import matplotlib.pyplot as plt
plt.style.use('bmh')
import argparse

# Currently only working in python 2 due to MDAnalysis package.
# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description = 'Run and score hdbscan clustering', add_help=False) 

# List all possible user input
inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-s', action='store', dest='structure',help='Structure file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-t', action='store', dest='trajectory',help='Trajectory',type=str,required=True)
inputs.add_argument('-sel', action='store', dest='sel', help='Atom selection',type=str,default='not element H')
inputs.add_argument('-min', action='store', dest='min', help='minimum cluster membership',type=int,default=10)
inputs.add_argument('-o', action='store', dest='out_name',help='Output file',type=str,required=True)

# Parse into useful form
UserInput=parser.parse_args()

topology = UserInput.structure
trajectory = UserInput.trajectory
t = md.load(trajectory,top=topology)
sel = t.topology.select(UserInput.sel)
t = t.atom_slice(sel)

out_name=UserInput.out_name
fname = './' + 'HD_' + out_name
os.mkdir(fname)

# Format trajectory 
temp = t.xyz
frames = t.xyz.shape[0]
atoms = t.xyz.shape[1]
data = temp.reshape((frames,atoms*3))
data = data.astype('float64')
temp = []

# Run hdbscan
clusterer = hdbscan.HDBSCAN(min_cluster_size=UserInput.min)
cluster_labels = clusterer.fit_predict(data)
raw_score = silhouette_score(data,cluster_labels)

# Predict QT cutoff
labels_reduced = cluster_labels[cluster_labels != -1]
counts = np.bincount(labels_reduced)
biggest_cluster = np.argmax(counts)
frames_big_cluster = np.where(cluster_labels == biggest_cluster)
biggest_cluster_traj = t.slice(frames_big_cluster)
# Make dissimilarity matrix for biggest cluster
distances = np.empty((biggest_cluster_traj.n_frames, biggest_cluster_traj.n_frames))
for i in range(biggest_cluster_traj.n_frames):
    distances[i] = md.rmsd(biggest_cluster_traj, biggest_cluster_traj, i)
print('Predicted QT diameter: %f nm' % np.max(distances))

# Calculate silhouette score ignoring points HDBSCAN treated as outliers
scoreable = np.where(cluster_labels>-1)[0]
scores = silhouette_samples(data,cluster_labels)
adjusted_score = np.mean(scores[scoreable])

# Save results
#Text
np.savetxt('HD_' + UserInput.out_name + '/hdbscan_labels.txt', cluster_labels, fmt='%i')
with open ('HD_' + UserInput.out_name + '/silhouette_scores.txt', 'w') as f:
    f.write("silhouette score is {0} \n adjusted score is {1}".format(raw_score, adjusted_score))
    
#Figures
plt.figure()
plt.scatter(np.arange(frames), cluster_labels, marker = '+')
plt.xlabel('Frame')
plt.ylabel('Cluster')
plt.title('HDBSCAN')
plt.savefig('HD_' + UserInput.out_name + '/hdbscan_timeseries.png')
plt.clf()
