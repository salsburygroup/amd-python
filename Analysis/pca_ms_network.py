import mdtraj as md
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import MeanShift
import networkx as nx
import matplotlib
matplotlib.use('Agg') #Avoids hangs on visualize
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description = 'Turn PCA space into a network based on wells', add_help=False) 

# List all possible user input
inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-top', action='store', dest='structure',help='Structure file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-traj', action='store', dest='trajectory',help='Trajectory',type=str,required=True)
inputs.add_argument('-sel', action='store', dest='sel', help='Atom selection',type=str,default='not element H')
inputs.add_argument('-o', action='store', dest='out_name',help='Output file',type=str,required=True)

# Parse into useful form
UserInput=parser.parse_args()

topology = UserInput.structure
trajectory = UserInput.trajectory
t = md.load(trajectory,top=topology)
sel = t.topology.select(UserInput.sel)
t = t.atom_slice(sel)

# Format trajectory 
temp = t.xyz
frames = t.xyz.shape[0]
atoms = t.xyz.shape[1]
data = temp.reshape((frames,atoms*3))
data = data.astype('float64')
temp = []

# Make PCA projection
projection = PCA().fit_transform(data)

num_components = projection.shape[1]
network = nx.Graph()
nodes = np.arange(0,num_components)
network.add_nodes_from(nodes)
label_dict = dict(zip(nodes,nodes))

for i in range(0,num_components-1):
    for j in range(i+1, num_components):
        temp = np.column_stack((projection[:,i], projection[:,j]))
        labels = MeanShift(n_jobs=-1).fit_predict(temp)
        num_clusters = max(labels) + 1 #0-based index
        if num_clusters > 1:
            network.add_edge(i,j)

nx.draw_spring(network, labels=label_dict)
plt.savefig(UserInput.out_name + '/PCA_network.png')
