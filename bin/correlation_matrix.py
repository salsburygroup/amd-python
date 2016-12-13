#!/usr/bin/env python

import argparse
from Analysis import AtomSelection, Correlation, Plotter, Saver, TrajectoryReader, TrajectoryProcessor


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
inputs.add_argument('-tau',
                    action='store',
                    dest='covariance_tau',
                    default=None,
                    type=int,
                    help='Lag time for constructing a time-lagged correlation matrix',
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

# Process trajectory
trajectory = TrajectoryReader.DCD(topology_path=UserInput.structure, trajectory_path=UserInput.trajectory).load()
trajectory = AtomSelection.Slice(trajectory=trajectory, atom_selection=UserInput.sel).select()

if UserInput.align:
    trajectory = TrajectoryProcessor.Aligner(trajectory=trajectory, atom_selection=UserInput.sel).process()

# Make correlation matrix

if UserInput.covariance_tau:
    correlation_matrix = Correlation.TimeLagged(
        trajectory=trajectory, covariance_tau=UserInput.covariance_tau
    ).calculate()
    title = 'Correlation Matrix with tau = {0}'.format(UserInput.covariance_tau)
else:
    correlation_matrix = Correlation.Pearson(trajectory=trajectory).calculate()
    title = 'Correlation Matrix'

# Save HeatMap

Plotter.UnityPColor(y=correlation_matrix,
                out_name=UserInput.out_name+'.png',
                x_label=UserInput.sel,
                y_label=UserInput.sel,
                title=title
                ).plot()

Saver.Array(
    array=correlation_matrix,
    out_name=UserInput.out_name+'.txt'
).save()
