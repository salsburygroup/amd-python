#! /usr/bin/env python

import argparse
import os
import numpy
import pickle
import mdtraj as md
import pandas as pd
from mdentropy.metrics import * 
from Analysis import Plotter, Saver

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
    description='Compute Solvent Accessible Surface Area  for given MD trajectories', add_help=False
)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-top',
                    action='store',
                    dest='structure',
                    help='Structure file corresponding to trajectory',
                    type=str,
                    required=True
                    )
inputs.add_argument('-traj',
                    action='store',
                    dest='trajectory',
                    help='Trajectory',
                    type=str,
                    required=True
                    )
inputs.add_argument('-stride',
                    action='store',
                    dest='stride',
                    help='stride to use',
                    type=int,
                    default=1
                    )
inputs.add_argument('-o',
                    action='store',
                    dest='out_name',
                    help='Output folder for text and png files',
                    type=str,
                    required=True
                    )
# Parse into useful form
UserInput = parser.parse_args()
traj = md.load(UserInput.trajectory, top=UserInput.structure, stride=UserInput.stride)
sasa = md.shrake_rupley(traj)
total_sasa = sasa.sum(axis=1)
                                   
Saver.Array(
    array=sasa,
    out_name=UserInput.out_name+'.txt'
).save()


title = 'Solvent Accessible Surface Area ({0})'.format(UserInput.out_name)
Plotter.Y(y=total_sasa,
                out_name=UserInput.out_name+'_total_sasa.pdf',
                x_label='Frames',
                y_label='Total SASA (nm)^2',
                title=title
                ).plot()

def autocorr(x):
    '''Compute an autocorrelation with numpy'''
    x = x - numpy.mean(x)
    result = numpy.correlate(x, x, mode='full')
    result = result[result.size//2:]
    return result / result[0]

autocorrelation = autocorr(total_sasa)
numpy.savetxt('autocorr.dat',autocorrelation)

import matplotlib
matplotlib.use('Agg') # For use on DEAC cluster
from matplotlib.pylab import *

semilogx(traj.time, autocorrelation)
xlabel('Time [ps]', size=16)
ylabel('SASA autocorrelation', size=16)
title('SASA autocorrelation ({0})'.format(UserInput.out_name))
savefig(UserInput.out_name+'_autocorr.pdf')
clf()

