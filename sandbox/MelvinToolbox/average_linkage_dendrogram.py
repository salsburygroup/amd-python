#!/usr/bin/env python

import argparse
import numpy
import matplotlib.pyplot as plt
import seaborn
import scipy.cluster.hierarchy
from scipy.cluster.hierarchy import dendrogram
import fastcluster
from scipy.spatial.distance import squareform

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
    description='Generate a publication-quality average linkage dendrogram from a distance matrix and ground truths.',
    add_help=False
)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-d',
                    action='store',
                    dest='distances',
                    help='Text file containing distance matrix',
                    type=str,
                    required=True
                    )
inputs.add_argument('-g',
                    action='store',
                    dest='ground_truths',
                    help='Text file containing ground-truth labels',
                    type=str,
                    default=None
                    )
inputs.add_argument('-t',
                    action='store',
                    dest='title',
                    help='Dendrogram title',
                    type=str,
                    default='Average linkage dendrogram'
                    )
inputs.add_argument('-x',
                    action='store',
                    dest='x_label',
                    help='Label for x-axis',
                    type=str,
                    default='sample index'
                    )
inputs.add_argument('-y',
                    action='store',
                    dest='y_label',
                    help='Label for y-axis',
                    type=str,
                    default='distance'
                    )
inputs.add_argument('-l',
                    action='store',
                    dest='labels',
                    help='Labels for x-ticks',
                    type=str,
                    default=None
                    )
inputs.add_argument('-c',
                    action='store',
                    dest='cutoff',
                    help='Where to cut the dendrogram',
                    type=float,
                    default=None
                    )
inputs.add_argument('-b',
                    action='store',
                    dest='bound',
                    help='Lower bound for y-axis',
                    type=float,
                    default=0.0
                    )
inputs.add_argument('-o',
                    action='store',
                    dest='out_name',
                    help='Output file name',
                    type=str,
                    required=True
                    )

# Parse into useful form
UserInput = parser.parse_args()

# Read in data
distance_matrix = numpy.genfromtxt(UserInput.distances)
distance_matrix = squareform(distance_matrix)
if UserInput.labels:
    with open(UserInput.labels) as l:
        #labels = l.read().strip().split(' ')
        labels = l.readlines()
else:
    labels = None


# Make linkage matrix
Z =fastcluster.linkage(distance_matrix, method='average')


scipy.cluster.hierarchy.set_link_color_palette(['orange', 'c', 'c', 'b', 'r'])


# Plot the dendrogram
fig = plt.figure()
d = dendrogram(Z, labels=labels, color_threshold=UserInput.cutoff)
leaves = d['leaves']
plt.title(UserInput.title)
plt.xlabel(UserInput.x_label)
plt.ylabel(UserInput.y_label)
ax = fig.gca()
ax.set_ylim(UserInput.bound, ax.get_ylim()[1])
ax.set_xticklabels(ax.xaxis.get_majorticklabels())

if UserInput.ground_truths:
    ground_truths = numpy.genfromtxt(UserInput.ground_truths).astype(int)
    #colors = plt.cm.rainbow(numpy.linspace(0, 1, max(ground_truths) + 1))
    colors = [(0, 0, 0), (.9, .6, 0), (.35, .7, .9), (0, .6, .5),
              (.95, .9, .25), (0, .46, .7), (.8, .4, 0), (.8, .6, .7)]
    ordered_ground_truths = ground_truths[leaves]
    idx = 0
    for x in ordered_ground_truths:
        plt.gca().get_xticklabels()[idx].set_color(colors[x])
        idx += 1
plt.savefig(UserInput.out_name, dpi=1200)
plt.close()
