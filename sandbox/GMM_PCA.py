import mdtraj as md
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import MeanShift
import networkx as nx
import matplotlib
matplotlib.use('Agg') #Avoids hangs on visualize
import matplotlib.pyplot as plt
import argparse
import sys

parser = argparse.ArgumentParser(description = 'Turn PCA space into a network based on wells', add_help=False) 

# List all possible user input
inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-top', action='store', dest='structure',help='Structure file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-traj', action='store', dest='trajectory',help='Trajectory',type=str,required=True)
inputs.add_argument('-sel', action='store', dest='sel', help='Atom selection',type=str,default='name CA')
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