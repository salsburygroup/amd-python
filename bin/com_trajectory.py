#!/usr/bin/env python

import argparse
import numpy
from Analysis import AtomSelection, Featurizer, Plotter, Saver, TrajectoryReader

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Calculate the center of mass for a selection in every frame',
                                 add_help=False
                                 )

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-s',
                    action='store',
                    dest='structure',
                    help='Structure file corresponding to trajectory',
                    type=str,
                    required=True
                    )
inputs.add_argument('-t',
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

com_trajectory = Featurizer.COM(trajectory=trajectory).extract()
Saver.Array(array=com_trajectory, out_name=UserInput.out_name + '_traj.txt').save()

com_change = [(com_frame - com_trajectory[0]) ** 2 for com_frame in com_trajectory]
com_change = numpy.sum(com_change, axis=-1)
com_change = numpy.sqrt(com_change)
Saver.Array(array=com_change, out_name=UserInput.out_name + '_chane.txt').save()
Plotter.Y(y=com_change,
          out_name=UserInput.out_name + '_change.png',
          x_label='Frame',
          y_label='COM change (nm)',
          title= 'COM Euclidean distance from starting frame').plot()
