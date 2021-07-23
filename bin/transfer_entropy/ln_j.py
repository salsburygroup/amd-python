#!/usr/bin/env python
from Analysis import TrajectoryReader, AtomSelection
from scipy.spatial.distance import cdist
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
p1=np.zeros((residues,bins))
for i in range(residues):
    p1[i]=np.histogram(delR[i], bins=bins)[0]/frames

# Calculating averaged ln of p1, p2, p2_tau, and p3
ln_j=np.nansum(p1*np.log(p1), axis=1)

# Saving p1, p2, p2_tau, and p3
np.save(os.path.join(UserInput.out_dir, 'ln_j.npy'), ln_j)
