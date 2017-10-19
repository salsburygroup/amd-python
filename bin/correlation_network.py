import numpy as np
import mdtraj as md
import argparse


# Jiajie Xiao
# 09.19.2017

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
    description='Compute temporal distribution of given time series of a variable', add_help=False
)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-input',
                    action='store',
                    dest='file',
                    help='input files of correlation matrix',
                    type=str,
                    required=True
                    )
inputs.add_argument('-str',
                    action='store',
                    dest='structure',
                    help='PDB file',
                    type=str,
                    required=True
                    )
inputs.add_argument('-o',
                    action='store',
                    dest='out_name',
                    help='Output file name',
                    type=str,
                    required=False
                    )

# Parse into useful form
UserInput = parser.parse_args()
topology = md.load(UserInput.structure).topology
residue_list = [residue for residue in topology.residues]

corr = np.genfromtxt(UserInput.file)
num_residues = corr.shape[0]

with open('corr_table.csv', 'w') as out:
    out.write('Res_A, Res_B, Type, Corr, abs(Corr)\n')
    for id_a in range(num_residues):
        for id_b in range(id_a+1, num_residues):
            out.write('{0:s}, {1:s}, {2:s}, {3:f}, {4:f}\n'.format(
                str(residue_list[id_a]), str(residue_list[id_b]), str(int(np.sign(corr[id_a, id_b]))),
                corr[id_a, id_b], abs(corr[id_a, id_b])))

