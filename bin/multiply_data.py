#!/usr/bin/env python
import argparse
import numpy as np
import pandas as pd
from Analysis import Saver

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Subtract two data', add_help=False)

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

inputs.add_argument('-m',
                    action='store',
                    dest='m',
                    help='multiply',
                    type=int,
                    default=1
                    )

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output name',
                    type=str,
                    default='multiply'
                    )

# Parse into useful form
UserInput = parser.parse_args()

data = pd.read_csv(UserInput.d,header=None)

pd.set_option('display.max_columns',999999,'display.max_rows',999999,)
print(data*UserInput.m)

# Save
Saver.Array(array=data*UserInput.m,
            out_name=UserInput.out_dir+'.dat'
            ).save()
