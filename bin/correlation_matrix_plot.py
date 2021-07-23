#!/usr/bin/env python
import os
import numpy as np
import argparse
from Analysis import Plotter

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Calculate, save and plot correlation matrix', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-d',
                    action='store',
                    dest='d',
                    help='data',
                    type=str,
                    required=True)

inputs.add_argument('-m',
                    action='store',
                    dest='m',
                    help='multiply',
                    type=int,
                    default=1)

inputs.add_argument('-title',
                    action='store',
                    dest='title',
                    help='Title of the plot',
                    type=str,
                    default=os.getcwd().split('/')[-1])

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output prefix for text and png',
                    type=str,
                    default='corr')

# Parse into useful form
UserInput = parser.parse_args()

# Make correlation matrix
correlation_matrix = np.loadtxt(UserInput.d)
title = 'Correlation Matrix (' + UserInput.title + ')'

# Plot HeatMap
Plotter.UnityPColor(y=correlation_matrix,
                    out_name=UserInput.out_dir,
                    x_label='residue number',
                    y_label='residue number',
                    title=title).plot()
