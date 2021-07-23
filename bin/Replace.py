#!/usr/bin/env python
import argparse

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Calculate, save and plot RMSF', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-fi',
                    action='store',
                    dest='fin',
                    help='Input file',
                    type=str,
                    required=True)

inputs.add_argument('-fo',
                    action='store',
                    dest='fout',
                    help='Output file',
                    type=str,
                    required=True)

inputs.add_argument('-ci',
                    action='store',
                    dest='cin',
                    help='Origin content',
                    type=str,
                    required=True)

inputs.add_argument('-co',
                    action='store',
                    dest='cout',
                    help='Replacing content',
                    type=str,
                    required=True)

# Parse into useful form
UserInput = parser.parse_args()

# Replacing contents
with open(UserInput.fin,'rt') as fin:
    with open(UserInput.fout,'wt') as fout:
        for line in fin:
            fout.write(line.replace(UserInput.cin, UserInput.cout))
