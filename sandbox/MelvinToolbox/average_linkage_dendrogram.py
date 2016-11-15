#!/usr/bin/env python

import argparse
import numpy
import matplotlib.pyplot as plt
import seaborn
from scipy.cluster.hierarchy import dendrogram, linkage

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
                    required=True
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
ground_truths = numpy.genfromtxt(UserInput.ground_truths).astype(int)

# Make linkage matrix
Z = linkage(distance_matrix, 'average')

# Plot the dendrogram
colors = plt.cm.rainbow(numpy.linspace(0,1,max(ground_truths)+1))
plt.figure()
# plt.style.use('ggplot')
plt.title(UserInput.title)
plt.xlabel(UserInput.x_label)
plt.xlabel(UserInput.y_label)
dendrogram(Z)
idx = 0
for x in ground_truths:
    plt.gca().get_xticklabels()[idx].set_color(colors[x])
    idx += 1
plt.savefig(UserInput.out_name, dpi=1200)
plt.close()
# plt.show()