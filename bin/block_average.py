import numpy as np
import matplotlib
matplotlib.use('Agg') # For use on DEAC cluster
import matplotlib.pyplot as plt
import argparse

# Jiajie Xiao
# 01.31.2017

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
    description='Compute standard errors of block average of given time series of a variable', add_help=False
)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-input',
                    action='store',
                    dest='file',
                    help='input files of a time series',
                    type=str,
                    required=True
                    )
inputs.add_argument('-min',
                    action='store',
                    dest='minBlockSize',
                    help='minimum block size',
                    type=int,
                    default=1
                    )
inputs.add_argument('-max',
                    action='store',
                    dest='maxBlockSize',
                    help='maximum block size',
                    type=int,
                    default=100
                    )
inputs.add_argument('-n',
                    action='store',
                    dest='numBlockSize',
                    help='number of block size to compute',
                    type=int,
                    default=50
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

timeSeries = np.genfromtxt(UserInput.file)

blockIncrement = 1 if 1.0*(UserInput.maxBlockSize - UserInput.minBlockSize)/UserInput.numBlockSize < 1 \
    else int((UserInput.maxBlockSize - UserInput.minBlockSize)/UserInput.numBlockSize)
blockSize_array = np.arange(UserInput.minBlockSize, UserInput.maxBlockSize, blockIncrement)

ste_list = []
for blockSize in blockSize_array:
    numBlock = int(len(timeSeries)/blockSize)
    blockAverage = np.mean(timeSeries[:int(numBlock*blockSize)].reshape(-1, int(blockSize)), 1)
    ste = blockAverage.std()/(numBlock**0.5)
    ste_list.append(ste)

np.savetxt(UserInput.out_name + '_bse.dat', np.column_stack([blockSize_array, ste_list]))

plt.figure()
plt.plot(blockSize_array, ste_list, '-b+')
plt.xlabel('Block Size')
plt.ylabel('BSE')
plt.title('Blocked Standard Error')
plt.savefig(UserInput.out_name + 'BSE.pdf')
plt.close()

