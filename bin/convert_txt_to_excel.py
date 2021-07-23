import pandas as pd
import argparse

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Calculate, save and plot correlation matrix', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-d',
                    action='store',
                    dest='d',
                    help='Input data',
                    type=str,
                    required=True
                    )

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output prefix for text and png',
                    type=str,
                    required=True
                    )

# Parse into useful form
UserInput = parser.parse_args()

df = pd.read_csv(UserInput.d, header=None, sep='\s+')
df.to_excel(UserInput.out_dir + '.xlsx', 'Sheet1')

