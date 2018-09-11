import numpy as np
import matplotlib
matplotlib.use('Agg') # For use on DEAC cluster
import matplotlib.pyplot as plt
import argparse

# Jiajie Xiao
# 01.31.2017

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
    description='Compute autocorrelation of given time series of a variable', add_help=False
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


def autocorr(x):
    "Compute an autocorrelation with numpy"
    x = x - np.mean(x)
    result = np.correlate(x, x, mode='full')
    result = result[result.size//2:]
    return result / result[0]

autoCorr = autocorr(timeSeries)
np.savetxt(UserInput.out_name + '_autocorrelation.dat', autoCorr)

#np.savetxt(UserInput.out_name + '.dat', np.column_stack([blockSize_array, ste_list]))

plt.figure()
plt.plot(range(len(autoCorr)), autoCorr, '-b+')
plt.xlabel('Time')
plt.ylabel('Autocorrelation')
plt.title('Autocorrelation')
plt.savefig(UserInput.out_name + '_autocorrelation.pdf')
plt.close()

