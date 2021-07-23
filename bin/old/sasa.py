#! /usr/bin/env python
import matplotlib as plt
plt.use('Agg') # For use on DEAC cluster
from matplotlib.pylab import *
from Analysis import Plotter, Saver
from mdentropy.metrics import * 
import mdtraj as md
import argparse
import numpy
import os

#Jiajie Xiao
# Jan 09, 2017
# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
    description='Compute Solvent Accessible Surface Area  for given MD trajectories', add_help=False
)

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

inputs.add_argument('-title',
                    action='store',
                    dest='title',
                    help='Title',
                    type=str,
                    required=True)

inputs.add_argument('-stride',
                    action='store',
                    dest='stride',
                    help='stride to use',
                    type=int,
                    default=1)

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output folder for text and png files',
                    type=str,
                    required=True)

# Parse into useful form
UserInput = parser.parse_args()
traj = md.load(UserInput.trajectory, top=UserInput.structure, stride=UserInput.stride)
sasa = md.shrake_rupley(traj)
total_sasa = sasa.sum(axis=1)
path = os.getcwd()

# Make output directory
fname = path + '/SASA'
if not os.path.exists(fname):
    os.mkdir(fname)                                

Saver.Array(
    array=sasa,
    out_name=os.path.join(fname, UserInput.out_dir + '.txt')
).save()

title = 'Solvent Accessible Surface Area ({0})'.format(UserInput.title)
Plotter.Y(y=total_sasa,
          out_name=os.path.join(fname, UserInput.out_dir + '_total_sasa.pdf'),
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
numpy.savetxt(os.path.join(fname, 'autocorr.dat'), autocorrelation)

fig2=plt.figure(2,figsize=(16,4))
semilogx(traj.time, autocorrelation)
plt.xlabel('Time [ps]', size=16)
plt.ylabel('SASA autocorrelation', size=16)
plt.title('SASA autocorrelation ({0})'.format(UserInput.title))
fig2.savefig(os.path.join(fname, UserInput.out_dir + '_autocorr.pdf'))
plt.clf()

