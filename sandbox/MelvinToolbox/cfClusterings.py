import numpy as np
import argparse
from sklearn import metrics

#Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description = 'Compare two clustering timeseries, with consideration for different labelings and number of clusters. 1 is best possible -1 is worst possible.', add_help=False) 

inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-t1', action='store', dest='timeseries1',help='First time series', type=str,required=True)
inputs.add_argument('-t2', action='store', dest='timeseries2', help='Second time series',required=True)

UserInput=parser.parse_args()

# Get time series 1
chain1 = np.genfromtxt(UserInput.timeseries1)

# Get time series 2
chain2 = np.genfromtxt(UserInput.timeseries2)

sm = metrics.adjusted_rand_score(chain1, chain2)

print('Adjusted rand index between two labelings is {}'.format(sm))

