#! /usr/bin/env python
import argparse
import numpy
import matplotlib.pyplot
import matplotlib.mlab
import seaborn

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
    description='Calculate RMSD time series. Centers structures at' +
                ' origin before beginning but does not otherwise align.',
    add_help=False
)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument(
    '-i',
    action='store',
    dest='input',
    help='series to histogram',
    type=str,
    required=True
)
inputs.add_argument(
    '-o',
    action='store',
    dest='output',
    help='output file name',
    type=str,
    required=True
)
inputs.add_argument(
    '-b',
    action='store',
    dest='bins',
    help='specify number, else FD rule will be used',
    type=int,
    default=None
)
# Parse into useful form
UserInput = parser.parse_args()

# Read input
timeseries = numpy.genfromtxt(fname=UserInput.input)
if not UserInput.bins:
    q75 = numpy.percentile(timeseries, 75)
    q25 = numpy.percentile(timeseries, 25)
    irq = q75 - q25
    bin_width = 2 * irq / (len(timeseries) ** (1. / 3))
    bins = int(numpy.ceil((numpy.max(timeseries) - numpy.min(timeseries)) / bin_width))
else:
    bins = UserInput.bins
n = matplotlib.pyplot.hist(x=timeseries, bins=bins, normed=1)
y = matplotlib.mlab.normpdf( bins, numpy.mean(timeseries), numpy.std(a=timeseries))
matplotlib.pyplot.grid(True)
matplotlib.pyplot.savefig(UserInput.output)
