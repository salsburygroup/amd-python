import numpy as np
import argparse

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Cropping images', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-d',
                    action='store',
                    dest='data',
                    help='Input data',
                    nargs='+',
                    required=True)

inputs.add_argument('-o',
                    action='store',
                    dest='out_name',
                    help='Output name',
                    type=str,
                    required=True)

# Parse into useful form
UserInput=parser.parse_args()

# Loading data
data=[]
for i in range(len(UserInput.data)):
    print(UserInput.data[i])
    data.append(np.loadtxt(UserInput.data[i]))
data=np.array(data)

# Calculating mean
mean=np.mean(data, axis=0)

# Saving mean
np.savetxt(UserInput.out_name, mean, delimiter=' ')

## Example
#python ~/python/Mean.py -d /deac/salsburyGrp/wud18/md/double_binding/1/transfer_entropy_max0.08/transfer_entropy_Na_stride1_1/transfer_entropy.dat /deac/salsburyGrp/wud18/md/double_binding/2/transfer_entropy_max0.08/transfer_entropy_Na_stride1_1/transfer_entropy.dat /deac/salsburyGrp/wud18/md/double_binding/3/transfer_entropy_max0.08/transfer_entropy_Na_stride1_1/transfer_entropy.dat /deac/salsburyGrp/wud18/md/double_binding/4/transfer_entropy_max0.08/transfer_entropy_Na_stride1_1/transfer_entropy.dat /deac/salsburyGrp/wud18/md/double_binding/5/transfer_entropy_max0.08/transfer_entropy_Na_stride1_1/transfer_entropy.dat /deac/salsburyGrp/wud18/md/double_binding/6/transfer_entropy_max0.08/transfer_entropy_Na_stride1_1/transfer_entropy.dat /deac/salsburyGrp/wud18/md/double_binding/7/transfer_entropy_max0.08/transfer_entropy_Na_stride1_1/transfer_entropy.dat /deac/salsburyGrp/wud18/md/double_binding/8/transfer_entropy_max0.08/transfer_entropy_Na_stride1_1/transfer_entropy.dat /deac/salsburyGrp/wud18/md/double_binding/9/transfer_entropy_max0.08/transfer_entropy_Na_stride1_1/transfer_entropy.dat /deac/salsburyGrp/wud18/md/double_binding/10/transfer_entropy_max0.08/transfer_entropy_Na_stride1_1/transfer_entropy.dat /deac/salsburyGrp/wud18/md/double_binding/11/transfer_entropy_max0.08/transfer_entropy_Na_stride1_1/transfer_entropy.dat /deac/salsburyGrp/wud18/md/double_binding/12/transfer_entropy_max0.08/transfer_entropy_Na_stride1_1/transfer_entropy.dat -o transfer_entropy_average1.dat
