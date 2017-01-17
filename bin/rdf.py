#! /usr/bin/env python

import argparse
import os
import numpy
import mdtraj as md
from Analysis import Plotter, Saver

# Jiajie Xiao
# 01.09.2017

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
    description='Compute radial distribution function of given selections for given MD trajectories', add_help=False
)

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
inputs.add_argument('-stride',
                    action='store',
                    dest='stride',
                    help='stride to use',
                    type=int,
                    default=1
                    )
inputs.add_argument('-minr',
                    action='store',
                    dest='minr',
                    help='minimum radius',
                    type=float,
                    default=0.0
                    )
inputs.add_argument('-maxr',
                    action='store',
                    dest='maxr',
                    help='maximum radius',
                    type=float,
                    default=1.0
                    )
inputs.add_argument('-periodic',
                    action='store',
                    dest='periodicTrue',
                    help='periodic true of false',
                    type=bool,
                    default=True
                    )
inputs.add_argument('-sel1',
                    action='store',
                    dest='sel1',
                    help='selection 1',
                    type=str,
                    required=True
                    )
inputs.add_argument('-sel2',
                    action='store',
                    dest='sel2',
                    help='selection 2',
                    type=str,
                    required=True
                    )
inputs.add_argument('-o',
                    action='store',
                    dest='out_name',
                    help='Output folder for text and png files',
                    type=str,
                    required=True
                    )
# Parse into useful form
UserInput = parser.parse_args()
traj = md.load(UserInput.trajectory, top=UserInput.structure, stride=UserInput.stride)
pais_Of_interest = traj.topology.select_pairs(selection1=UserInput.sel1,selection2=UserInput.sel2)
radii_range = (UserInput.minr, UserInput.maxr)

r, g_r = md.compute_rdf(traj,pairs=pais_Of_interest,r_range=radii_range,bin_width=0.01,n_bins=None,periodic=UserInput.periodicTrue,opt=True)
rdf = numpy.column_stack((r,g_r))                                  
Saver.Array(
    array=rdf,
    out_name='rdf_'+UserInput.out_name+'.dat'
).save()


title = 'Radial distribution function ({0},{1})'.format(UserInput.sel1, UserInput.sel2)
Plotter.XY(x=r,y=g_r,
                out_name='rdf_'+UserInput.out_name+'.pdf',
                x_label='r (nm)',
                y_label='g(r)',
                title=title
                ).plot()
