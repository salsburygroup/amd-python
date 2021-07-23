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

# Calculating net entropy transfer from each residue
E_net=np.sum(Tij-Tij.T, axis=1)

# Saving transfer entropy and net entropy transfer 
np.savetxt(os.path.join(UserInput.out_dir, UserInput.nm + '.dat'), Tij, delimiter=' ')
np.savetxt(os.path.join(UserInput.out_dir, UserInput.nm + '_net.dat'), E_net, delimiter=' ')

# Plot HeatMap
if UserInput.max == None:
    vmax=np.max(Tij)
else:
    vmax=UserInput.max
    
fig1=plt.figure(1,figsize=(8,6))
plt.pcolor(Tij.T, cmap='jet', vmin=0, vmax=vmax)
cbar=plt.colorbar()
plt.title(UserInput.title)
plt.xlabel(UserInput.axis_X)
plt.ylabel(UserInput.axis_Y)
fig1.savefig(os.path.join(UserInput.out_dir, UserInput.nm + '.png'), dpi=200)
fig1.savefig(os.path.join(UserInput.out_dir, UserInput.nm + '.tiff'), dpi=600)
fig1.savefig(os.path.join(UserInput.out_dir, UserInput.nm + '.pdf'), bbox_inches='tight')

fig2=plt.figure(2,figsize=(16,4))
plt.plot(E_net)
plt.title(UserInput.title)
plt.xlabel('residue number')
plt.ylabel('net entropy transfer')
fig2.savefig(os.path.join(UserInput.out_dir, UserInput.nm + '_net.png'), dpi=200)
fig2.savefig(os.path.join(UserInput.out_dir, UserInput.nm + '_net.tiff'), dpi=600)
fig2.savefig(os.path.join(UserInput.out_dir, UserInput.nm + '_net.pdf'), bbox_inches='tight')
