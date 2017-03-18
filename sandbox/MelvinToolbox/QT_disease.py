import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg') #For use on DEAC cluster
import matplotlib.pyplot as plt
plt.style.use('bmh')
import tempfile
import argparse
from scipy.spatial.distance import squareform

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description = 'Runs a vectorized version of QT clustering', add_help=False) 

# List all possible user input
inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-min', action='store', dest='minimum_membership', help='Minimum number of frames in a cluster',type=int,default=2)
inputs.add_argument('-cutoff', action='store', dest='cutoff', help='maximum cluster radius',type=float,required=True)
inputs.add_argument('-distance', action='store', dest='distances', help='distance matrix',type=str,required=True)
inputs.add_argument('-labels', action='store', dest='labels', help='labels for points',type=str, default=None)
inputs.add_argument('-o', action='store', dest='out_name',help='Output directory',type=str,required=True)

# Parse into useful form
UserInput = parser.parse_args()
distances = np.genfromtxt(UserInput.distances)
if UserInput.labels:
    with open(UserInput.labels) as l:
        x_labels = l.read().strip().split(' ')
else:
    x_labels = None


cutoff_mask = distances <= UserInput.cutoff
centers = []
cluster = 0
labels = np.empty(distances.shape[0])
labels.fill(np.NAN)

while cutoff_mask.any():
    membership = cutoff_mask.sum(axis=1)
    center = np.argmax(membership)
    members = np.where(cutoff_mask[center,:]==True)
    if max(membership) <= UserInput.minimum_membership:
        labels[np.where(np.isnan(labels))] = -1
        break
    labels[members] = cluster
    centers.append(center)
    cutoff_mask[members,:] = False
    cutoff_mask[:,members] = False
    cluster = cluster + 1


# Save results
#Text
np.savetxt(UserInput.out_name + '/QT_labels.txt', labels, fmt='%i')
np.savetxt(UserInput.out_name + '/QT_centers.txt', centers, fmt='%i')

#Figures
fig = plt.figure()
plt.scatter(np.arange(distances.shape[0]), labels, marker = '+')
ax = fig.gca()
ax.set_xticklabels(x_labels)
plt.xlabel('Disease ID')
plt.ylabel('Cluster')
plt.title('QT')
plt.savefig(UserInput.out_name + '/QT.png')
plt.close()
