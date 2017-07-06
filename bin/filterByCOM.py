#! /usr/bin/env python
import argparse
import numpy
from copy import deepcopy
import sys
import mdtraj
from Analysis import Distance, Plotter, Saver, TrajectoryReader, TrajectoryProcessor, Featurizer, AtomSelection

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
    description='Calculate COM time series, between two separate atom selections.' +
                'Useful for determining if a ligand dissociates from binding pocket,' +
                'for example',
    add_help=False
)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument(
    '-structure',
    action='store',
    dest='structure',
    help='Structure file corresponding to trajectory',
    type=str,
    required=True
)
inputs.add_argument('-stride',
                    action='store',
                    dest='stride',
                    help='stride to use',
                    type=int,
                    default=1
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
    '-sel1',
    action='store',
    dest='sel1',
    help='Atom selection',
    type=str,
    default='all'
)

inputs.add_argument(
    '-sel2',
    action='store',
    dest='sel2',
    help='Atom selection 2',
    type=str,
    default='all'
)
inputs.add_argument(
    '-cutoff_low',
    action='store',
    dest='cutoffl',
    help='Distance cutoff (In nanometers)',
    type=float,
    default=5.0
)
inputs.add_argument(
    '-cutoff_high',
    action='store',
    dest='cutoffh',
    help='Distance cutoff (In nanometers)',
    type=float,
    default=5.0
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

# Read trajectory
trajectory_all = mdtraj.load(UserInput.trajectory, top=UserInput.structure, stride=UserInput.stride)
#trajectory_all = TrajectoryReader.DCD(topology_path=UserInput.structure, trajectory_path=UserInput.trajectory).load()

#print(UserInput.sel1)

sub_trajectory1 = AtomSelection.Slice(trajectory=trajectory_all, atom_selection=UserInput.sel1).select()
sub_trajectory2 = AtomSelection.Slice(trajectory=trajectory_all, atom_selection=UserInput.sel2).select()

# Calculate Center Of Mass of each atom selection and the difference thereof
COM1=Featurizer.COM(sub_trajectory1).extract()
COM2=Featurizer.COM(sub_trajectory2).extract()

dists = (COM1 - COM2)**2
distances = dists.sum(axis=-1)
COM_diff=numpy.sqrt(distances)
print(COM_diff)

#Create a mask where the cutoff condition is satisfied
masky_low = numpy.where(COM_diff<UserInput.cutoffl)
masky_med = numpy.where((COM_diff>=UserInput.cutoffl) & (COM_diff<=UserInput.cutoffh))
masky_high = numpy.where(COM_diff>UserInput.cutoffh)

print(masky_high)

masked_traj_low = trajectory_all[masky_low]
print(masked_traj_low)

masked_traj_med = trajectory_all[masky_med]
print(masked_traj_med)

masked_traj_high = trajectory_all[masky_high]
print(masked_traj_high)
#print(masked_COM)

# Save trajectories
masked_traj_low.save_dcd(filename=UserInput.out_name + '_low.dcd')
masked_traj_med.save_dcd(filename=UserInput.out_name + '_med.dcd')
masked_traj_high.save_dcd(filename=UserInput.out_name + '_high.dcd')

# Plot
Plotter.Y(
    y=COM_diff,
    out_name=UserInput.out_name + '.png',
    x_label='Frame',
    y_label='COM Diff (nm)',
    title='COM Distance timeseries'
).plot()
