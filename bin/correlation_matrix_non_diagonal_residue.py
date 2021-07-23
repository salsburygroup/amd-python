#!/usr/bin/env python
import argparse
import mdtraj as md
import numpy as np

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Calculate, save and plot correlation matrix', add_help=False)

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

topology = md.load(UserInput.structure).topology
residue_list = [residue for residue in topology.residues]
corr = np.genfromtxt(UserInput.in_dir + '.dat')
num_residues = corr.shape[0]

with open((UserInput.out_dir + '_non_diagnonal.dat'), 'w') as out:
    for id_a in range(num_residues):
        id_b_start = id_a + 1 
        for id_b in range(id_b_start, num_residues):
            out.write('{0:f}\n'.format(corr[id_a, id_b]))

with open((UserInput.out_dir + '_non_diagnonal_residue.dat'), 'w') as out:
    out.write('Res_A, Res_B, Type, Corr\n')
    for id_a in range(num_residues):
        id_b_start = id_a + 1 
        for id_b in range(id_b_start, num_residues):
            out.write('{0:s}, {1:s}, {2:s}, {3:f}\n'.format(str(residue_list[id_a]), str(residue_list[id_b]), str(int(np.sign(corr[id_a, id_b]))), corr[id_a, id_b]))
