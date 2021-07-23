#!/usr/bin/env python
import argparse
import numpy as np

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Calculate, save and plot correlation matrix', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-i',
                    action='store',
                    dest='in_dir',
                    help='Input prefix for text and png',
                    type=str,
                    default='corr'
                    )

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output prefix for text and png',
                    type=str,
                    default='corr'
                    )

# Parse into useful form
UserInput = parser.parse_args()

corr = np.genfromtxt(UserInput.in_dir + '.dat')
num_residues = corr.shape[0]

with open((UserInput.out_dir + '_non_diagnonal.dat'), 'w') as out:
    for id_a in range(num_residues):
        id_b_start = id_a + 1
        for id_b in range(id_b_start, num_residues):
            out.write('{0:f}\n'.format(corr[id_a, id_b]))
