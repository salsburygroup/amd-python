#!/usr/bin/env python

import argparse
from Analysis.TrajectoryProcessor import Repair


parser = argparse.ArgumentParser(
    description='Repair a dcd that has atoms within residues in a different order than intended',
    add_help=False
                                 )

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-w',
                    action='store',
                    dest='wrong_structure',
                    help='Structure file with incorrect atom ordering (but same residue ordering)',
                    type=str,
                    required=True
                    )
inputs.add_argument('-r',
                    action='store',
                    dest='right_structure',
                    help='Structure file with correct atom ordering',
                    type=str,
                    required=True
                    )
inputs.add_argument('-t',
                    action='store',
                    dest='trajectory_file',
                    help='Trajectory file corresponding to wrong_structure',
                    type=str,
                    required=True
                    )
inputs.add_argument('-o',
                    action='store',
                    dest='out_name',
                    help='Output name for repaired DCD',
                    type=str,
                    required=True
                    )

# Parse into useful form
UserInput = parser.parse_args()

Repair.atom_order(
    right_pdb_file=UserInput.right_structure,
    wrong_pdb_file=UserInput.wrong_structure,
    wrong_trajectory_file=UserInput.trajectory_file,
    out_trajectory_file=UserInput.out_name
)
