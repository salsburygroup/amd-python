#!/usr/bin/env python
from Analysis import TrajectoryReader, AtomSelection
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import argparse
import os

print('Make sure to use the absolute path for output directory(-o)')
# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Compute the transfer entropy for the selected residues', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-title',
                    action='store',
                    dest='title',
                    help='Title of the plot',
                    type=str,
                    default='transfer entropy')

inputs.add_argument('-axX',
                    action='store',
                    dest='axis_X',
                    default='residue number (entropy donor)',
                    help='Label for X axes')

inputs.add_argument('-axY',
                    action='store',
                    dest='axis_Y',
                    default='residue number (entropy acceptor)',
                    help='Label for Y axes')

inputs.add_argument('-nm',
                    action='store',
                    dest='nm',
                    help='Output name for data',
                    type=str,
                    default='transfer_entropy')

inputs.add_argument('-max',
                    action='store',
                    dest='max',
                    help='Upper limit for the HeatMap',
                    type=float,
                    default=None)

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output folder for data',
                    type=str,
                    default=os.getcwd())

# Parse into useful form
UserInput = parser.parse_args()

# Loading the averaged ln of p1, p2, p2_tau, and p3
ln_jjtau=np.load('ln_jjtau.npy')
ln_ijjtau=np.load('ln_ijjtau.npy')
ln_j=np.load('ln_j.npy')
ln_ij=np.load('ln_ij.npy')

# Calculating transfer entropy
Tij=-ln_jjtau+ln_ijjtau+ln_j-ln_ij
        
# Saving transfer entropy
np.savetxt(os.path.join(UserInput.out_dir, UserInput.nm + '.dat'), Tij, delimiter=' ')

# Plot HeatMap
if UserInput.max == True:
    vmax=UserInput.max
else:
    vmax=np.max(Tij)
    
fig1=plt.figure(1,figsize=(8,6))
plt.pcolor(Tij.T, cmap='jet', vmin=0, vmax=vmax)
cbar=plt.colorbar()
plt.title(UserInput.title)
plt.xlabel(UserInput.axis_X)
plt.ylabel(UserInput.axis_Y)
fig1.savefig(os.path.join(UserInput.out_dir, UserInput.nm + '.png'), dpi=200)
fig1.savefig(os.path.join(UserInput.out_dir, UserInput.nm + '.tiff'), dpi=600)
fig1.savefig(os.path.join(UserInput.out_dir, UserInput.nm + '.pdf'), bbox_inches='tight')

fig2=plt.figure(2,figsize=(8,4))
plt.plot(-ln_jjtau)
plt.title(UserInput.title)
plt.xlabel(UserInput.axis_X)
plt.ylabel(UserInput.axis_Y)
fig2.savefig(os.path.join('ln_jjtau.png'), dpi=200)

fig3=plt.figure(3,figsize=(8,6))
plt.pcolor(ln_ijjtau, cmap='jet', vmin=np.min(ln_ijjtau), vmax=np.max(ln_ijjtau))
cbar=plt.colorbar()
plt.title(UserInput.title)
plt.xlabel(UserInput.axis_X)
plt.ylabel(UserInput.axis_Y)
fig3.savefig(os.path.join('ln_ijjtau.png'), dpi=200)

fig4=plt.figure(4,figsize=(8,4))
plt.plot(ln_j)
plt.title(UserInput.title)
plt.xlabel(UserInput.axis_X)
plt.ylabel(UserInput.axis_Y)
fig4.savefig(os.path.join('ln_j.png'), dpi=200)

fig5=plt.figure(5,figsize=(8,6))
plt.pcolor(-ln_ij, cmap='jet', vmin=np.min(-ln_ij), vmax=np.max(-ln_ij))
cbar=plt.colorbar()
plt.title(UserInput.title)
plt.xlabel(UserInput.axis_X)
plt.ylabel(UserInput.axis_Y)
fig5.savefig(os.path.join('ln_ij.png'), dpi=200)


