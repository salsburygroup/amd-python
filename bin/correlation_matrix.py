#!/usr/bin/env python

import argparse
from Analysis import AtomSelection, Correlation, Plotter, Saver, TrajectoryReader

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Calculate, save and plot correlation matrix', add_help=False)

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
inputs.add_argument('-o',
                    action='store',
                    dest='out_name',
                    help='Output prefix for text and png',
                    type=str,
                    required=True
                    )

# Parse into useful form
UserInput = parser.parse_args()

# Process trajectory
trajectory = TrajectoryReader.DCD(topology_path=UserInput.structure, trajectory_path=UserInput.trajectory).load()
trajectory = AtomSelection.Slice(trajectory=trajectory, atom_selection=UserInput.sel).select()

# Make correlation matrix
correlation_matrix = Correlation.Pearson(trajectory=trajectory).calculate()

# Save HeatMap
Plotter.SimplePColor(y=correlation_matrix,
                out_name=UserInput.out_name+'.png',
                x_label=UserInput.sel,
                y_label=UserInput.sel,
                title='Correlation Matrix'
                ).plot()

Saver.Array(
    array=correlation_matrix,
    out_name=UserInput.out_name+'.txt'
).save()
