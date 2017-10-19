#! /usr/bin/env python

import argparse
import os
import numpy as np
import mdtraj as md
import pandas as pd
import matplotlib
matplotlib.use('Agg') # For use on DEAC cluster
import matplotlib.pyplot as plt
#from Analysis import Plotter, Cluster


# Jiajie Xiao
# 07.19.2017

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
    description='Compute temporal distribution of given time series of a variable', add_help=False
)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-input',
                    action='store',
                    nargs='+',
                    dest='files',
                    help='input files containing the time series of a variable',
                    type=str,
                    required=True
                    )
inputs.add_argument('-binWidth',
                    action='store',
                    dest='bin_width',
                    help='width of the bins for variable histogramming and discretizing ',
                    type=float,
                    required=True
                    )
inputs.add_argument('-o',
                    action='store',
                    dest='out_name',
                    help='Output folder for data',
                    type=str,
                    required=True
                    )

# Parse into useful form
UserInput = parser.parse_args()

# Load input files to the variable x
x = np.array([np.loadtxt(file) for file in UserInput.files])
num_files, num_samples = x.shape

# Histogramming and discretize
left = x.min()
right = x.max()
num_bins = int((right - left)/UserInput.bin_width)
bin_edges = np.linspace(left, right, num_bins+1)
inds = np.digitize(x, bin_edges) # bin indices of each data point, starting from 1 in our case
hist, _ = np.histogram(x,bins=bin_edges,normed=False) # the numpy histogram with a normalization is buggy

plt.figure()
plt.plot(.5*(bin_edges[1:]+bin_edges[:-1]), hist,'-*')
plt.xlabel('Variable')
plt.ylabel('Count')
plt.title('Histogram of ({0})'.format(UserInput.out_name))
plt.savefig('hist_'+UserInput.out_name+'bin='+str(UserInput.bin_width)+'.pdf')
plt.close()

plt.figure()
plt.scatter(np.arange(inds.size), inds.ravel(), marker='+')
plt.xlabel('Frames')
plt.ylabel('Bin id')
plt.title('Timeseries of discretized {0}'.format(UserInput.out_name))
plt.savefig('TimeSeries_'+UserInput.out_name+'bin='+str(UserInput.bin_width)+'.pdf')
plt.close()

#Plotter.XY(x=.5*(bin_edges[1:]+bin_edges[:-1]), y=hist,
#                out_name='hist_'+UserInput.out_name+'bin='+str(UserInput.bin_width)+'.pdf',
#                x_label='Variable',
#                y_label='Probability',
#                title='Histogram of ({0})'.format(UserInput.out_name)
#                ).plot()
#Cluster.Plotter.TimeSeries(inds,'TimeSeries_'+UserInput.out_name+'bin='
#                                         +str(UserInput.bin_width)+'.pdf').plot()

# Residence and off time
#continuity = np.zeros(x.shape, dtype=int) # matrix representing if the current value is same as the previous value; 1 denotes yes; 0 denotes no
residence_time = np.zeros(x.shape, dtype=int)
for run in range(num_files):
    for index in range(num_samples-1):
        if inds[run,-(index+1)] == inds[run,-index]: # within the same bin
            #continuity[run,-(index+1)]=1
            residence_time[run,-(index+1)] = residence_time[run,-index]+1
        else:
            residence_time[run,-(index+1)] = 1 # the minimum residence time should be 1
    residence_time[run, num_samples-1] = 1 # assume a transition would occur after the last observation, this may result in some errors
    residence_time[run, 0] = 1
# mean residence time
mean_residence_time = [residence_time.ravel()[(ind+1)==inds.ravel()].mean() for ind in range(num_bins+1)]

#Plotter.XY(x=bin_edges, y=mean_residence_time,
#                out_name='TemporalDistribution'+UserInput.out_name+'bin='+str(UserInput.bin_width)+'.pdf',
 #               x_label='Variable',
 #               y_label='Frames',
 #               title='TemporalDistribution of ({0})'.format(UserInput.out_name)
 #               ).plot()

plt.figure()
plt.plot(bin_edges, mean_residence_time,'-*')
plt.xlabel('Variable')
plt.ylabel('Number of Frames')
plt.title('TemporalDistribution of {0}'.format(UserInput.out_name))
plt.savefig('TemporalDistribution'+UserInput.out_name+'bin='+str(UserInput.bin_width)+'.pdf')
plt.close()

temporal_distribution = pd.DataFrame(np.column_stack((bin_edges, mean_residence_time)), columns=['Variable values','Residence time'])
temporal_distribution.to_csv('TemporalDistribution'+UserInput.out_name+'bin='+str(UserInput.bin_width)+'.csv',header=True, index=False)







