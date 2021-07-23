from Analysis.Cluster import Saver
import argparse
import glob
import os

print('Make sure to use the absolute path for output directory(-o)')

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Visualizing clusters', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-o', action='store', dest='out_dir', help='Absolute path of output directory', type=str, default=os.getcwd())

# Parse into useful form
UserInput = parser.parse_args()

# Visualizing clusters
for folder in glob.glob(os.path.join(UserInput.out_dir, 'clusters', 'cluster*')):
    Saver.Shadows(out_name=folder, middle=os.path.join(folder, 'rep.pdb'), shadow=os.path.join(folder, 'all.dcd')).save()
    #Saver.Shadows(out_name=folder, middle=os.path.join(folder, 'rep.pdb'), shadow=os.path.join(folder, 'all.pdb')).save()
