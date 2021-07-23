#!/usr/bin/env python
from Analysis import TrajectoryReader, AtomSelection
from scipy.spatial.distance import cdist
import mdtraj as md
import numpy as np
import subprocess
import argparse
import os

print('Make sure to use the absolute path for output directory(-o)')
# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Compute the transfer entropy for the selected residues', add_help=False)

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

inputs.add_argument('-sel',
                    action='store',
                    dest='sel',
                    help='Atom selection',
                    type=str,
                    default='name CA')

inputs.add_argument('-title',
                    action='store',
                    dest='title',
                    help='Title of the plot',
                    type=str,
                    default='transfer\\\ entropy')

inputs.add_argument('-l',
                    action='store',
                    dest='stride',
                    help='Stride to use',
                    type=str,
                    required=True)

inputs.add_argument('-tau',
                    action='store',
                    dest='tau',
                    help='Decay time (ns)',
                    type=str,
                    default=1)

inputs.add_argument('-axX',
                    action='store',
                    dest='axis_X',
                    default='residue\\\\\ number\\\\\ \\\\\(entropy\\\\\ donor\\\\\)',
                    help='Label for X axes')

inputs.add_argument('-axY',
                    action='store',
                    dest='axis_Y',
                    default='residue\\\\\ number\\\\\ \\\\\(entropy\\\\\ acceptor\\\\\)',
                    help='Label for Y axes')

inputs.add_argument('-nm',
                    action='store',
                    dest='nm',
                    help='Output name for data',
                    type=str,
                    default='transfer_entropy')

inputs.add_argument('-max',
                    action='store',
                    dest='max',
                    help='Upper limit for the HeatMap',
                    type=str,
                    default='None')

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output folder for data',
                    type=str,
                    default='transfer_entropy')

# Parse into useful form
UserInput = parser.parse_args()

# Make output directory
if not os.path.exists(UserInput.out_dir):
    os.mkdir(UserInput.out_dir)

# Calculating the averaged value of ln separately
subprocess.call([os.getenv('SHELL'), '-i', '-c', 'cd ' + UserInput.out_dir + ' && python3 /home/wud18/python/transfer_entropy/ln_ijTrial.py -s ' + UserInput.structure + ' -t ' + UserInput.trajectory + ' -sel ' + UserInput.sel + ' -l ' + UserInput.stride + ' -tau ' + UserInput.tau + ' -o ' + UserInput.out_dir + ' & sleep 0.01 && python3 /home/wud18/python/transfer_entropy/ln_jjtauTrial.py -s ' + UserInput.structure + ' -t ' + UserInput.trajectory + ' -sel ' + UserInput.sel + ' -l ' + UserInput.stride + ' -tau ' + UserInput.tau + ' -o ' + UserInput.out_dir + ' & sleep 0.01 && python3 /home/wud18/python/transfer_entropy/ln_jTrial.py -s ' + UserInput.structure + ' -t ' + UserInput.trajectory + ' -sel ' + UserInput.sel + ' -l ' + UserInput.stride + ' -tau ' + UserInput.tau + ' -o ' + UserInput.out_dir + ' & sleep 0.1 && python3 /home/wud18/python/transfer_entropy/ln_ijjtauTrial.py -s ' + UserInput.structure + ' -t ' + UserInput.trajectory + ' -sel ' + UserInput.sel + ' -l ' + UserInput.stride + ' -tau ' + UserInput.tau + ' -title ' + UserInput.title + ' -axX ' + UserInput.axis_X + ' -axY ' + UserInput.axis_Y + ' -nm ' + UserInput.nm + ' -max ' + UserInput.max + ' -o ' + UserInput.out_dir + ' &'], cwd=UserInput.out_dir)

## Using example

#python3 ~/python/transfer_entropy/transfer_entropy.py -s /deac/salsburyGrp/wud18/md/double_binding/protein_Na_right_angle.pdb -t /deac/salsburyGrp/wud18/md/double_binding/protein_Na_stride100_aligned.dcd -sel name\\\\\\\ CA\\\\\\\ or\\\\\\\ resname\\\\\\\ SOD -l 100 -tau 6 -o /deac/salsburyGrp/wud18/md/double_binding/transfer_entropy_Na_stride100_6 -title transfer\\\\\\\ entropy\\\\\\\ \\\\\\\(tau=6\\\\\\\)
