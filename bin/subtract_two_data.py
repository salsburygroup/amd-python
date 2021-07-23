#!/usr/bin/env python
from Analysis import Saver
import pandas as pd
import argparse

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Subtract two data', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-d1',
                    action='store',
                    dest='d1',
                    help='data1',
                    type=str,
                    required=True)

inputs.add_argument('-d2',
                    action='store',
                    dest='d2',
                    help='data2',
                    type=str,
                    required=True)

inputs.add_argument('-m',
                    action='store',
                    dest='m',
                    help='multiply',
                    type=int,
                    default=1)

inputs.add_argument('-o',
                    action='store',
                    dest='out_name',
                    help='Output name',
                    type=str,
                    default='subtract')

# Parse into useful form
UserInput = parser.parse_args()

data1 = pd.read_csv(UserInput.d1,header=None, sep='\s+')
data2 = pd.read_csv(UserInput.d2,header=None, sep='\s+')

pd.set_option('display.max_columns',999999,'display.max_rows',999999,)
print((data1 - data2)*UserInput.m)

# Save
Saver.Array(array=(data1-data2)*UserInput.m,
            out_name=UserInput.out_name+'.dat'
            ).save()
