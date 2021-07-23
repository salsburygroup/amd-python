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
inputs.add_argument('-s',
                    action='store',
                    dest='structure',
                    help='Structure file corresponding to trajectory',
                    type=str,
                    required=True)

inputs.add_argument('-t',
                    action='store',
                    dest='trajectory',
                    help='Trajectory',
                    type=str,
                    required=True)

inputs.add_argument('-sel',
                    action='store',
                    dest='sel',
                    help='Atom selection',
                    type=str,
                    default='name CA')

inputs.add_argument('-title',
                    action='store',
                    dest='title',
                    help='Title of the plot',
                    type=str,
                    default='transfer entropy')

inputs.add_argument('-l',
                    action='store',
                    dest='stride',
                    help='Stride to use',
                    type=int,
                    required=True)

inputs.add_argument('-tau',
                    action='store',
                    dest='tau',
                    help='Decay time (ns)',
                    type=int,
                    default=1)

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
                    default='transfer_entropy')

# Parse into useful form
UserInput = parser.parse_args()

# Make output directory
if not os.path.exists(UserInput.out_dir):
    os.mkdir(UserInput.out_dir)

# Loading trajectory
trajectory = TrajectoryReader.DCD(topology_path=UserInput.structure, trajectory_path=UserInput.trajectory).load()

# Calculating the averaged structures
sub_trajectory = AtomSelection.Slice(trajectory=trajectory, atom_selection=UserInput.sel).select()
assert isinstance(sub_trajectory, md.Trajectory)
averaged_structures = sub_trajectory.xyz.mean(0)

# Calculating the number of bins
frames=len(sub_trajectory)
residues=len(averaged_structures)
bins=int(1+np.log(frames)/np.log(2))
tau=int(100*UserInput.tau/UserInput.stride)
steps=frames-tau

# Calculating fluctuations
tmp1=(sub_trajectory.xyz-averaged_structures).reshape(frames*residues,3)
tmp2=cdist(tmp1,[[0,0,0]],metric='euclidean')
delR=tmp2.reshape(frames,residues)
delR=delR.T

# Binning p1, p2, p2_tau, and p3
p1=np.zeros((residues,bins))
for i in range(residues):
    p1[i]=np.histogram(delR[i], bins=bins)[0]/frames

p2=np.zeros((residues,residues,bins,bins))
for i in range(residues):
    for j in range(residues):
        p2[i][j]=np.histogramdd([delR[i], delR[j]], bins=bins)[0]/frames

p2_tau=np.zeros((residues,residues,bins,bins))
if tau == 0 :
    for i in range(residues):
        for j in range(residues):
            p2_tau[i][j]=np.histogramdd([delR[i], delR[j]], bins=bins)[0]/steps
else:
    for i in range(residues):
        for j in range(residues):
            p2_tau[i][j]=np.histogramdd([delR[i][:-tau], delR[j][tau:]], bins=bins)[0]/steps

p3=np.zeros((residues,residues,bins,bins,bins))
if tau == 0 :
    for i in range(residues):
        for j in range(residues):
            p3[i][j]=np.histogramdd([delR[i], delR[j], delR[j]], bins=bins)[0]/steps
else:
    for i in range(residues):
        for j in range(residues):
            p3[i][j]=np.histogramdd([delR[i][:-tau], delR[j][:-tau], delR[j][tau:]], bins=bins)[0]/steps

# Calculating averaged ln of p1, p2, p2_tau, and p3
ln_jjtau=np.nansum(np.nansum(p2_tau*np.log(p2_tau), axis=2), axis=2)
ln_ijjtau=np.nansum(np.nansum(np.nansum(p3*np.log(p3), axis=2), axis=2), axis=2)
ln_j=np.nansum(p1*np.log(p1), axis=1)
ln_ij=np.nansum(np.nansum(p2*np.log(p2), axis=2), axis=2)

# Saving p1, p2, p2_tau, and p3
np.save(os.path.join(UserInput.out_dir, 'ln_jjtau.npy'), ln_jjtau)
np.save(os.path.join(UserInput.out_dir, 'ln_ijjtau.npy'), ln_ijjtau)
np.save(os.path.join(UserInput.out_dir, 'ln_j.npy'), ln_j)
np.save(os.path.join(UserInput.out_dir, 'ln_ij.npy'), ln_ij)

# Calculating transfer entropy
Tij=-ln_jjtau.T+ln_ijjtau.T+ln_j-ln_ij.T
        
# Saving transfer entropy
np.savetxt(os.path.join(UserInput.out_dir, UserInput.nm + '.dat'), Tij.T, delimiter=' ')

# Plot HeatMap
if UserInput.max == True:
    vmax=UserInput.max
else:
    vmax=np.max(Tij)
    
fig1=plt.figure(1,figsize=(8,6))
plt.pcolor(Tij, cmap='jet', vmin=0, vmax=vmax)
cbar=plt.colorbar()
plt.title(UserInput.title)
plt.xlabel(UserInput.axis_X)
plt.ylabel(UserInput.axis_Y)
fig1.savefig(os.path.join(UserInput.out_dir, UserInput.nm + '.png'), dpi=200)
fig1.savefig(os.path.join(UserInput.out_dir, UserInput.nm + '.tiff'), dpi=600)
fig1.savefig(os.path.join(UserInput.out_dir, UserInput.nm + '.pdf'), bbox_inches='tight')

# Using example
#python3 ~/python/transfer_entropy_combined.py -s protein_Na_right_angle.pdb -t protein_Na_stride100_aligned.dcd -sel 'name CA or resname SOD' -o transfer_entropy_Na_stride100_2 -tau 2 -l 100
