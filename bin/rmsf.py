#!/usr/bin/env python

import argparse
from Analysis import Distance, Plotter, Saver, TrajectoryReader, TrajectoryProcessor

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Calculate, save and plot RMSF', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-top',
                    action='store',
                    dest='structure',
                    help='Structure file corresponding to trajectory',
                    type=str,
                    required=True
                    )
inputs.add_argument('-traj',
                    action='store',
                    dest='trajectory',
                    help='Trajectory',
                    type=str,
                    required=True
                    )
inputs.add_argument('-sel',
                    action='store',
                    dest='sel',
                    help='Atom selection',
                    type=str,
                    default='name CA'
                    )
inputs.add_argument('-align',
                    action='store_true',
                    help='Align to atom selection before calculating?',
)
inputs.add_argument('-o',
                    action='store',
                    dest='out_name',
                    help='Output prefix for text and png',
                    type=str,
                    required=True
                    )

# Parse into useful form
UserInput = parser.parse_args()

# Load trajectory
trajectory = TrajectoryReader.DCD(topology_path=UserInput.structure,
                                  trajectory_path=UserInput.trajectory
                                  ).load()

if UserInput.align:
    trajectory = TrajectoryProcessor.Aligner(trajectory=trajectory, atom_selection=UserInput.sel).process()

# Calculate RMSF
rmsf_array = Distance.RMSF(trajectory=trajectory,
                           atom_selection=UserInput.sel
                           ).calculate()

# Plot
Plotter.Y(y=rmsf_array,
          out_name=UserInput.out_name+'.png',
          x_label=UserInput.sel,
          y_label='RMSF (nm)',
          title='RMSF'
          ).plot()

# Save
Saver.Array(array=rmsf_array,
            out_name=UserInput.out_name+'.txt'
            ).save()
