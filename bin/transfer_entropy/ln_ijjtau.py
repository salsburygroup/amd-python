#!/usr/bin/env python
from Analysis import TrajectoryReader, AtomSelection
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import argparse
import time
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

inputs.add_argument('-l',
                    action='store',
                    dest='stride',
                    help='Stride to use',
                    type=str,
                    required=True)

inputs.add_argument('-tau',
                    action='store',
                    dest='tau',
                    help='Decay time (ns)',
                    type=str,
                    default=1)
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
                    type=str,
                    default='None')

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output folder for data',
                    type=str,
                    default='transfer_entropy')

# Parse into useful form
UserInput = parser.parse_args()

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
tau=int(100*float(UserInput.tau)/int(UserInput.stride))
steps=frames-tau

# Calculating fluctuations
tmp1=(sub_trajectory.xyz-averaged_structures).reshape(frames*residues,3)
tmp2=cdist(tmp1,[[0,0,0]],metric='euclidean')
delR=tmp2.reshape(frames,residues)
delR=delR.T

# Binning p1, p2, p2_tau, and p3
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
ln_ijjtau=np.nansum(np.nansum(np.nansum(p3*np.log(p3), axis=2), axis=2), axis=2)
time.sleep(1800)
ln_jjtau=np.load(os.path.join(UserInput.out_dir, 'ln_jjtau.npy'))
ln_j=np.load(os.path.join(UserInput.out_dir, 'ln_j.npy'))
ln_ij=np.load(os.path.join(UserInput.out_dir, 'ln_ij.npy'))

# Saving p1, p2, p2_tau, and p3
np.save(os.path.join(UserInput.out_dir, 'ln_ijjtau.npy'), ln_ijjtau)
np.savetxt(os.path.join(UserInput.out_dir, 'ln_ijjtau.txt'), ln_ijjtau)
np.savetxt(os.path.join(UserInput.out_dir, 'ln_jjtau.txt'), ln_jjtau)
np.savetxt(os.path.join(UserInput.out_dir, 'ln_j.txt'), ln_j)
np.savetxt(os.path.join(UserInput.out_dir, 'ln_ij.txt'), ln_ij)

# Calculating transfer entropy
Tij=-ln_jjtau+ln_ijjtau+ln_j-ln_ij

# Calculating net entropy transfer from each residue
E_net=np.sum(Tij-Tij.T, axis=1)

# Saving transfer entropy and net entropy transfer 
np.savetxt(os.path.join(UserInput.out_dir, UserInput.nm + '.dat'), Tij, delimiter=' ')
np.savetxt(os.path.join(UserInput.out_dir, UserInput.nm + '_net.dat'), E_net, delimiter=' ')

# Plot HeatMap
if UserInput.max == 'None':
    vmax=np.max(Tij)
else:
    vmax=UserInput.max
    
fig1=plt.figure(1,figsize=(8,6))
plt.pcolor(Tij.T, cmap='jet', vmin=0, vmax=vmax)
cbar=plt.colorbar()
plt.title(UserInput.title)
plt.xlabel(UserInput.axis_X)
plt.ylabel(UserInput.axis_Y)
fig1.savefig(os.path.join(UserInput.out_dir, UserInput.nm + '.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
fig1.savefig(os.path.join(UserInput.out_dir, UserInput.nm + '.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
fig1.savefig(os.path.join(UserInput.out_dir, UserInput.nm + '.pdf'), bbox_inches='tight')

fig2=plt.figure(2,figsize=(16,4))
plt.plot(E_net)
plt.title(UserInput.title)
plt.xlabel('residue number')
plt.ylabel('net entropy transfer')
fig2.savefig(os.path.join(UserInput.out_dir, UserInput.nm + '_net.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
fig2.savefig(os.path.join(UserInput.out_dir, UserInput.nm + '_net.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
fig2.savefig(os.path.join(UserInput.out_dir, UserInput.nm + '_net.pdf'), bbox_inches='tight')

# Plot each ln
fig3=plt.figure(3,figsize=(8,6))
plt.pcolor(ln_ijjtau, cmap='jet', vmin=np.min(ln_ijjtau), vmax=np.max(ln_ijjtau))
cbar=plt.colorbar()
plt.title('ln_ijjtau')
plt.xlabel(UserInput.axis_X)
plt.ylabel(UserInput.axis_Y)
fig3.savefig(os.path.join(UserInput.out_dir, 'ln_ijjtau.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)

fig4=plt.figure(4,figsize=(16,4))
plt.plot(ln_j)
plt.title('ln_j')
plt.xlabel('residue number')
plt.ylabel('ln_j')
fig4.savefig(os.path.join(UserInput.out_dir, 'ln_j.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)

fig5=plt.figure(5,figsize=(8,6))
plt.pcolor(ln_ij, cmap='jet', vmin=np.min(ln_ij), vmax=np.max(ln_ij))
cbar=plt.colorbar()
plt.title('ln_ij')
plt.xlabel(UserInput.axis_X)
plt.ylabel(UserInput.axis_Y)
fig5.savefig(os.path.join(UserInput.out_dir, 'ln_ij.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)

fig6=plt.figure(6,figsize=(16,4))
plt.plot(ln_jjtau)
plt.title('ln_jjtau')
plt.xlabel('residue number')
plt.ylabel('ln_jjtau')
fig6.savefig(os.path.join(UserInput.out_dir, 'ln_jjtau.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
