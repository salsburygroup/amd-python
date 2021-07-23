#!/usr/bin/env python
import os
import numpy as np
import argparse

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
                    required=True
                    )

inputs.add_argument('-th',
                    action='store',
                    dest='th',
                    help='Threshold for correlation coefficients',
                    type=float,
                    default=0.8
                    )

inputs.add_argument('-ig',
                    action='store',
                    dest='ig',
                    help='Ignore adjacent residues',
                    type=int,
                    default=1
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

corr = np.genfromtxt(UserInput.d)
num_residues = corr.shape[0]

# Correlation coefficients map to structure
with open(UserInput.out_dir + '_map.dat', 'w') as out:
    for id_a in range(num_residues):
        id_b_start = id_a
        for id_b in range(id_b_start, num_residues):
            out.write('{0:s}\t{1:s}\t{2:f}\n'.format(str(id_a), str(id_b), corr[id_a, id_b]))

# Correlation coefficients map to structure(corr > UserInput.th)
with open(UserInput.out_dir + '_map_' + str(UserInput.th) + '_ig' + str(UserInput.ig) + '.dat', 'w') as out:
    for id_a in range(num_residues):
        id_b_start = id_a + UserInput.ig
        for id_b in range(id_b_start, num_residues):
            if (abs(corr[id_a, id_b]) > UserInput.th) and (id_b>185 and id_a<186):
                out.write('{0:s}\t{1:s}\t{2:f}\n'.format(str(id_a), str(id_b), corr[id_a, id_b]))
