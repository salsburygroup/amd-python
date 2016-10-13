#! /usr/bin/env python
#
#psfgenWrapper.py
#Ryan Godwin
#9/30/16

import argparse
from Analysis import psfgen, Saver, TrajectoryReader, TrajectoryProcessor

parser = argparse.ArgumentParser(
    description='Reads a trajectory, makes a pdb for every frame'+
                ', apply psfgen to each pdb, making new psf '+
                ' and pdb, and combines new pdbs to new traj',
    epilog='Psfgen wrapping complete',
    add_help=False
)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument(
    '-struct',
    action='store',
    dest='structure',
    help='Structure file corresponding to trajectory',
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
    '-top',
    action='store',
    dest='topology',
    help='Topology',
    type=str,
    required=True
)
inputs.add_argument(
    '-patch',
    action='store',
    nargs = 3,
    dest='patch',
    help='Patch',
    type=str,
    required=False
)
inputs.add_argument(
    '-o',
    action='store',
    dest='out_name',
    help='Output prefix',
    type=str,
    required=True
)

# Parse into useful form
UserInput = parser.parse_args()

# Read trajectory,  deteremine frame count, and apply seletcion mask
trajectory = TrajectoryReader.DCD(topology_path=UserInput.structure, trajectory_path=UserInput.trajectory).load()
n_frames = trajectory.n_frames
#sel = t.topology.select(UserInput.sel)
#t = t.atom_slice(sel)

#Check for Patches
if UserInput.patch:
    patch_type = UserInput.patch[0]
    patch_segment = UserInput.patch[1]
    patch_resid = int(UserInput.patch[2])
    print(patch_type, patch_segment, patch_resid)
else:
    print('No patch defined...')


print(UserInput.topology)


#topo_test = Psfgen(topology_path = UserInput.topology).readCharmmTopology()

#for loop over number of frames
for frame in range(n_frames):
    print(frame)
    #save pdb of frame
    traj_frame = trajectory.slice(frame)
    
    
    #generate a psf and pdb from that frame with correct parameters
    topo_new=psfgen(topology_path)

        
    

#psfgen
#recombine pdb files 
