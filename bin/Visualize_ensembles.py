from Analysis.Cluster import Saver_common
import argparse
import glob
import os

print('Make sure to use the absolute path for output directory(-o)')
# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Visualizing ensembles', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-s', action='store', dest='structure', help='Structure file corresponding to trajectory', type=str, required=True)
inputs.add_argument('-t', action='store', dest='trajectory', help='Trajectory', type=str, required=True)
inputs.add_argument('-o', action='store', dest='out_dir', help='Absolute path of output directory', type=str, required=True)

# Parse into useful form
UserInput = parser.parse_args()

# Make output directory
out_dir=UserInput.out_dir
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# Visualizing clusters
Saver_common.Shadows(out_name=UserInput.out_dir, middle=UserInput.structure, shadow=UserInput.trajectory).save()

