import pyemma.coordinates as coor
import numpy as np
import argparse

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Calculate, save and plot RMSF', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-s',
                    action='store',
                    dest='structure',
                    help='Structure file corresponding to trajectory',
                    type=str,
                    required=True)

inputs.add_argument('-t',
                    action='store',
                    dest='trajectory',
                    help='Trajectory',
                    type=str,
                    required=True)

inputs.add_argument('-n',
                    action='store',
                    dest='num',
                    help='Number of random frames',
                    type=str,
                    default='50')

inputs.add_argument('-o',
                    action='store',
                    dest='out_name',
                    help='Output name for .dcd',
                    type=str,
                    default='random.dcd')

# Parse into useful form
UserInput = parser.parse_args()

# Loading structure and trajectory
topfile = UserInput.structure
trajfile = UserInput.trajectory
feat = coor.featurizer(topfile)
inp = coor.source(trajfile, feat)
print('trajectory length = ',inp.trajectory_length(0))
print('number of dimension = ',inp.dimension())

# Selection 50 random frames
num=int(UserInput.num)
random=np.ceil(np.random.random(num)*num)-1
random_process=[]
for i in range(len(random)):
    random_process.append([0,int(random[i])])

# Saving random frames
coor.save_traj(inp, random_process, UserInput.out_name)
