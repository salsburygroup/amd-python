#! /usr/bin/env python
from Analysis import Distance, Saver, TrajectoryReader, TrajectoryProcessor
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

# Initialize parser.
parser = argparse.ArgumentParser(
    description='Plot rmsd and rmsf for all simulations', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-title',
                    action='store',
                    dest='title',
                    help='Title of the plot',
                    type=str,
                    default=os.getcwd().split('/')[-1])

inputs.add_argument('-n',
                    action='store',
                    dest='number',
                    help='The number of runs',
                    type=int,
                    default=5)

inputs.add_argument('-sel',
                    action='store',
                    dest='sel',
                    help='Atom selection',
                    type=str,
                    default='name CA')

inputs.add_argument('-tm',
                    action='store',
                    dest='timestep',
                    help='Simulation timestep (ps)',
                    type=int,
                    default=4)

inputs.add_argument('-fq',
                    action='store',
                    dest='dcdfreq',
                    help='dcd frequency',
                    type=int,
                    default='2500')

inputs.add_argument('-dpi',
                    action='store',
                    dest='dpi',
                    help='Dots per inch (resolution of the picture)',
                    type=int,
                    default=200)

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output prefix',
                    type=str,
                    default='CA')

# Parse into useful form
UserInput = parser.parse_args()

# Make output directory
out_dir=UserInput.out_dir
fname = './' + 'deviation_' + out_dir
if not os.path.exists(fname):
    os.mkdir(fname)

# Plot rmsd
t = np.arange(0, len(np.loadtxt('1/deviation_CA/rmsd_CA.dat'))*UserInput.dcdfreq*UserInput.timestep/1000000, UserInput.dcdfreq*UserInput.timestep/1000000)
fig1=plt.figure(1,figsize=(16,4))
for n in range(UserInput.number):
    plt.plot(t, np.loadtxt(str(n+1)+'/deviation_CA/rmsd_CA.dat')*10, alpha=0.9, linewidth=0.5, label='run'+str(n+1))

plt.legend(loc='lower right')
plt.xlabel('time (ns)')
plt.ylabel('RMSD (Å)')
plt.title(UserInput.title + ' RMSD timeseries')
fig1.savefig(os.path.join('deviation_'+out_dir,'rmsd_'+out_dir+'.png'), pad_inches=0.03, bbox_inches='tight', dpi=UserInput.dpi)
fig1.savefig(os.path.join('deviation_'+out_dir,'rmsd_'+out_dir+'.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
fig1.savefig(os.path.join('deviation_'+out_dir,'rmsd_'+out_dir+'.pdf'), bbox_inches='tight')

# Xlabel
if UserInput.sel == 'name CA':
    xlabel = 'residue number'
elif UserInput.sel == 'protein': 
    xlabel = 'atom number'
else:
    xlabel = UserInput.sel

# Plot rmsf
fig2=plt.figure(2,figsize=(16,4))
for n in range(UserInput.number):
    plt.plot(np.loadtxt(str(n+1)+'/deviation_CA/rmsf_CA.dat')*10, label='run'+str(n+1), alpha=0.9, linewidth=0.5)

plt.legend(loc='upper right')
plt.xlabel(xlabel)
plt.ylabel('RMSF (Å)')
plt.title('RMSF of ' + UserInput.title)
fig2.savefig(os.path.join('deviation_'+out_dir,'rmsf_'+out_dir+'.png'), pad_inches=0.03, bbox_inches='tight', dpi=UserInput.dpi)
fig2.savefig(os.path.join('deviation_'+out_dir,'rmsf_'+out_dir+'.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
fig2.savefig(os.path.join('deviation_'+out_dir,'rmsf_'+out_dir+'.pdf'), bbox_inches='tight')
