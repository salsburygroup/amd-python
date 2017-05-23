#! /usr/bin/env python
import argparse
from Analysis import Distance, Plotter, Saver, TrajectoryReader, TrajectoryProcessor

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
    description='Calculate RMSD time series. Centers structures at' +
                ' origin before beginning but does not otherwise align.',
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
inputs.add_argument(
    '-ref',
    action='store',
    dest='ref_frame',
    help='Reference Structure',
    type=int,
    default=0
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
    help='Atom selection',
    type=str,
    default='all'
)
inputs.add_argument(
    '-align',
    action='store_true',
    help='Align to atom selection before calculating?',
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
trajectory = TrajectoryReader.DCD(topology_path=UserInput.structure, trajectory_path=UserInput.trajectory).load()

if UserInput.align:
    trajectory = TrajectoryProcessor.Aligner(trajectory=trajectory, atom_selection=UserInput.sel).process()

# Calculate RMSD
rmsd_timeseries = Distance.RMSD(
    trajectory=trajectory, atom_selection=UserInput.sel, reference_frame=UserInput.ref_frame
).calculate()

# Save
Saver.Array(array=rmsd_timeseries, out_name=UserInput.out_name + '.dat').save()

# Plot
Plotter.Y(
    y=rmsd_timeseries,
    out_name=UserInput.out_name + '.png',
    x_label='Frame',
    y_label='RMSD (nm)',
    title='RMSD timeseries'
).plot()
