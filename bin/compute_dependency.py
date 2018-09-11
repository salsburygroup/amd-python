#! /usr/bin/env python
import argparse
import pyemma.coordinates as coor
import mdtraj as md
import numpy as np
from sklearn.feature_selection import mutual_info_regression
import matplotlib
matplotlib.use('Agg') # For use on DEAC cluster
import seaborn as sns; sns.set()
#%matplotlib inline


# Jiajie Xiao
# April 3 2018

# Example call
# python compute_dependency.py -str 4dih_fill.B99990001_autopsf.pdb -tr protein_all_aligned_stride10.dcd -sel "not element H" -o test_mutant_average

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
    description='Compute the residue-residue dependency based on the center of mass of each residue',
    add_help=False
)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument(
    '-str',
    action='store',
    dest='structure',
    help='PDB file not psf',
    type=str,
    required=True
)
inputs.add_argument(
    '-traj',
    action='store',
    dest='trajectory',
    help='Trajectory',
    type=str,
    required=True
)
inputs.add_argument(
    '-sel',
    action='store',
    dest='sel',
    help='Atom selection (Not implemented yet)',
    type=str,
    default='all'
)
inputs.add_argument(
    '-stride',
    action='store',
    dest='stride',
    help='striding',
    type=int,
    default='1'
)
inputs.add_argument(
    '-m',
    action='store',
    dest='method',
    help='pearson, covariance, mi or nmi',
    type=str,
    default='pearson'
)
inputs.add_argument(
    '-o',
    action='store',
    dest='out_name',
    help='Output file name',
    type=str,
    required=True
)

# Parse into useful form
UserInput = parser.parse_args()

traj = UserInput.trajectory
top = UserInput.structure
t = md.load(top)

num_residues = t.topology.n_residues
feat = coor.featurizer(top)
feat.add_residue_COM(range(num_residues))

X = coor.load(traj, feat,stride=UserInput.stride)
num_frames = X.shape[0]
Y = X.reshape([num_frames,num_residues,3]) # Tensor form
Y -= Y.mean(axis=0) # 0-means

def compute_dependency(x1, x2, method):
    if method == 'pearson' or method == 'covariance':
        return (x1*x2).sum(axis=1).mean()
    elif method == 'mi' or method == 'nmi':
        #return mutual_info_regression(Y[:,res_i,0].reshape(-1, 1),Y[:,res_j,1])+mutual_info_regression(Y[:,res_i,0].reshape(-1, 1),Y[:,res_j,2])+mutual_info_regression(Y[:,res_i,1].reshape(-1, 1),Y[:,res_j,2])
        mi = 0
        for n in range(x2.shape[1]):
            mi += (mutual_info_regression(Y[:,res_i,:],Y[:,res_j,n])).sum()
        return mi

dependency = np.zeros((num_residues,num_residues))
normalized_dependency = np.zeros((num_residues,num_residues))

for res_i in range(num_residues):
    for res_j in range(res_i,num_residues):
        dependency[res_i,res_j] = compute_dependency(Y[:,res_i,:], Y[:,res_j,:], UserInput.method)
        dependency[res_j,res_i] = dependency[res_i,res_j]

for res_i in range(num_residues):
    for res_j in range(res_i,num_residues):
        normalized_dependency[res_i,res_j] = dependency[res_i,res_j]/((dependency[res_i,res_i]*dependency[res_j,res_j])**0.5)
        normalized_dependency[res_j,res_i] = normalized_dependency[res_i,res_j]

# Extract maximum and minimum of the non-diagonal elements
mask = np.ones(dependency.shape, dtype=bool)
np.fill_diagonal(mask, 0)

nvMin = normalized_dependency[mask].min()
nvMin = normalized_dependency[mask].max()
vMin = dependency[mask].min()
vMax = dependency[mask].max()

fig = matplotlib.pyplot.figure()
#if UserInput.method == 'pearson':
#    vMin = -1;
#    vMax = 1;
#elif UserInput.method  == 'nmi':   
#    #vMin = 0;
#    #vMax = 1;
#    vMin = normalized_dependency[mask].min()
#    vMax = normalized_dependency[mask].max()
#elif UserInput.method  ==  ('mi' or 'covariance'):
#    #vMin = None
#    #vMax = None 
#    vMin = dependency[mask].min()
#    vMax = dependency[mask].max()


ax = sns.heatmap(dependency,xticklabels=50,yticklabels=50, 
                 vmin=vMin, vmax=vMax, cmap="YlGnBu",
                 #cbar_kws={'label': 'Success Rate of Climbs (%)'},
                 square=True)
ax.invert_yaxis()
ax.set_xlabel('Residue (COM)')
ax.set_ylabel('Residue (COM)')
if UserInput.method == 'pearson' or UserInput.method == 'covariance':
    ax.set_title('Cov ('+UserInput.out_name + ') ')
    fig.savefig('Cov_'+UserInput.out_name +'.pdf')
    np.savetxt('Cov_'+UserInput.out_name +'.txt', dependency, delimiter=',') 
elif UserInput.method  == 'mi' or UserInput.method  == 'nmi':
    ax.set_title('MI ('+UserInput.out_name + ') ')
    fig.savefig('MI_'+UserInput.out_name +'.pdf')
    np.savetxt('Mi_'+UserInput.out_name +'.txt', dependency, delimiter=',') 

fig = matplotlib.pyplot.figure()
ax = sns.heatmap(normalized_dependency,xticklabels=50,yticklabels=50, 
                 vmin=nvMin, vmax=nvMax, cmap="YlGnBu",
                 #cbar_kws={'label': 'Success Rate of Climbs (%)'},
                 square=True)
ax.invert_yaxis()
ax.set_xlabel('Residue (COM)')
ax.set_ylabel('Residue (COM)')
if UserInput.method == 'pearson' or UserInput.method == 'covariance':
    ax.set_title('Corr ('+UserInput.out_name + ') ')
    fig.savefig('Corr_'+UserInput.out_name +'.pdf')
    np.savetxt('Corr_'+UserInput.out_name +'.txt', normalized_dependency, delimiter=',') 
elif UserInput.method  == 'mi' or UserInput.method == 'nmi':
    ax.set_title('NMI ('+UserInput.out_name + ') ')
    fig.savefig('NMI_'+UserInput.out_name +'.pdf')
    np.savetxt('NMI_'+UserInput.out_name +'.txt', normalized_dependency, delimiter=',') 



