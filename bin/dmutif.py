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
    description='Compute mutual information for given MD trajectories', add_help=False
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
inputs.add_argument('-n',
                    action='store',
                    dest='nbins',
                    help='Number of bins',
                    type=int,
                    default=3
                    )
inputs.add_argument('-m',
                    action='store',
                    dest='method',
                    help='Entropy estimate method',
                    choices=['chaowangjost', 'grassberger', 'kde',
                            'knn', 'naive'],
                    default='knn'
                    )                    
inputs.add_argument('-it',
                    action='store',
                    dest='iter',
                    help='Number of shuffle iterations',
                    type=int,
                    default=100
                    )                    
inputs.add_argument('-threads',
                    action='store',
                    dest='n_threads',
                    help='Number of threads to be used',
                    type=int,
                    default=None
                    )
inputs.add_argument('-normed',
                    action='store',
                    dest='normed',
                    help='normalized mutual information or mutual information',
                    type=str,
                    default=False
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
mi = DihedralMutualInformation(n_bins=UserInput.nbins, method=UserInput.method, threads=UserInput.n_threads,
                                   normed=UserInput.normed)
                                   
M = mi.partial_transform(traj, shuffle=UserInput.iter, verbose=True)
df = pd.DataFrame(M, columns=mi.labels)
pickle.dump(df, open(UserInput.out_name, 'wb'))

title = 'Mutual information ({0})'.format(UserInput.method)
Plotter.UnityPColor(y=df,
                out_name=UserInput.out_name+'.png',
                x_label='Residue number',
                y_label='Residue number',
                title=title
                ).plot()

Saver.Array(
    array=df,
    out_name=UserInput.out_name+'.txt'
).save()
