#! /usr/bin/env python
import matplotlib
matplotlib.use('Agg') # For use on DEAC cluster
import matplotlib.pyplot as plt
from Analysis import Saver
import mdtraj as md
import argparse
import numpy
import os

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
    description='Compute Solvent Accessible Surface Area for given MD trajectories', add_help=False
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
                    required=True)

inputs.add_argument('-sel',
                    action='store',
                    dest='sel',
                    help='selected residues',
                    nargs='+',
                    default=None)

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output folder for text and png files',
                    type=str,
                    required=True)

# Parse into useful form
UserInput = parser.parse_args()

# Make output directory
out_dir = UserInput.out_dir
fname = out_dir
if not os.path.exists(fname):
    os.mkdir(fname)

# Calculate SASA
traj = md.load(UserInput.trajectory, top=UserInput.structure)
sasa = md.shrake_rupley(traj, mode='residue')
sasa_T = numpy.transpose(sasa)

total_sasa = sasa_T[int(UserInput.sel[0])]
for i in range(1,len(UserInput.sel)):
    total_sasa = total_sasa + sasa_T[int(UserInput.sel[i])]

total_sasa_mean=numpy.mean(total_sasa)
total_sasa_std=numpy.std(total_sasa)

numpy.savetxt(os.path.join(UserInput.out_dir, 'Total_sasa.dat'), total_sasa)
numpy.savetxt(os.path.join(UserInput.out_dir, 'Total_sasa_mean&std.dat'), [total_sasa_mean,total_sasa_std])

# Plot SASA
fig1=plt.figure(1,figsize=(16,4))
plt.plot(traj.time*UserInput.stride*0.01, total_sasa)
plt.xlabel('time(ns)')
plt.ylabel('Total SASA (nm)^2')
plt.title('Solvent Accessible Surface Area ({0})'.format(UserInput.title))
fig1.savefig(os.path.join(UserInput.out_dir, 'total_sasa.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
fig1.savefig(os.path.join(UserInput.out_dir, 'total_sasa.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
fig1.savefig(os.path.join(UserInput.out_dir, 'total_sasa.pdf'), bbox_inches='tight')
plt.clf()

# Calculate autocorr
def autocorr(x):
    '''Compute an autocorrelation with numpy'''
    x = x - numpy.mean(x)
    result = numpy.correlate(x, x, mode='full')
    result = result[result.size//2:]
    return result / result[0]

autocorrelation = autocorr(total_sasa)
numpy.savetxt(os.path.join(UserInput.out_dir, 'Autocorr.dat'), autocorrelation)

# Plot autocorr
fig2=plt.figure(2,figsize=(16,4))
plt.semilogx(traj.time*UserInput.stride*0.01, autocorrelation)
plt.xlabel('time(ns)')
plt.ylabel('SASA autocorrelation')
plt.title('SASA autocorrelation ({0})'.format(UserInput.title))
fig2.savefig(os.path.join(UserInput.out_dir, 'autocorr.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
fig2.savefig(os.path.join(UserInput.out_dir, 'autocorr.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
fig2.savefig(os.path.join(UserInput.out_dir, 'autocorr.pdf'), bbox_inches='tight')
