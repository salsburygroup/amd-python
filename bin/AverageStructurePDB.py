#! /usr/bin/env python
import argparse
import mdtraj
import numpy
from Analysis import TrajectoryReader, TrajectoryProcessor, Distance


# Jiajie Xiao
# Jun 18 2017

# Example call
# python AverageStructurePDB.py -str 4dih_fill.B99990001_autopsf.pdb -tr protein_all_aligned_stride10.dcd -sel "not element H" -o test_mutant_average

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
    description='Out put a PDB that is closest to the average structure (of interest) accross a trajectory',
    add_help=False
)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument(
    '-s',
    action='store',
    dest='structure',
    help='PDB file',
    type=str,
    required=True
)
inputs.add_argument(
    '-t',
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
    help='Atom selection',
    type=str,
    default='all'
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

# Read trajectory
#trajectory = TrajectoryReader.BigDCD(
#    topology_path=UserInput.structure,
#    trajectory_path=UserInput.trajectory,
#    chunk_size=100,
#).load()

#indices = AtomIndexer(self.trajectory, self.sel).process()
traj_temp = mdtraj.load_pdb(UserInput.structure)
indices = traj_temp.topology.select(UserInput.sel)
xyz_mean = traj_temp.xyz[0]
# Save the PDB with the averaged structure
xyz_mean_sel = numpy.zeros((len(indices),3))
num_frames = 0
for chunk in mdtraj.iterload(UserInput.trajectory,top=UserInput.structure,chunk=100):
    xyz_mean_sel += chunk.atom_slice(indices).xyz.sum(axis=0)
    num_frames += chunk.n_frames
xyz_mean_sel /= num_frames
xyz_mean[indices] = xyz_mean_sel

with mdtraj.formats.PDBTrajectoryFile(UserInput.out_name+'_average.pdb', 'w') as output_file:
    output_file.write(10*xyz_mean, traj_temp.topology)

# Calculate RMSD and find the frame that is closest average structure
rmsd_timeseries = []
traj_ref = mdtraj.load_pdb(UserInput.out_name + '_average.pdb')
# for chunk in trajectory:
for chunk in mdtraj.iterload(UserInput.trajectory,top=UserInput.structure, chunk=100):
    # rmsd_timeseries_temp = Distance.RMSD(
    #    trajectory=chunk, atom_selection=UserInput.sel, reference_frame=UserInput.out_name + '_average.pdb'
    # ).calculate()
    rmsd_timeseries.append(mdtraj.rmsd(chunk, traj_ref, frame=0, atom_indices=indices))
    # print(len(rmsd_timeseries))
average_frame_index = numpy.argmin(numpy.concatenate(rmsd_timeseries))
average_structure = mdtraj.load_frame(UserInput.trajectory, average_frame_index, top=UserInput.structure)
average_structure.save_pdb(UserInput.out_name+'_closest_average.pdb', 'w')
