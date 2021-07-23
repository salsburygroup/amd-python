import argparse
import os
import numpy
from Analysis import AtomSelection, DimensionReduction, Featurizer, FreeEnergy, Plotter, Saver, TrajectoryReader

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Calculate, save and plot principal components and projections', add_help=False)

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

inputs.add_argument('-n',
                    action='store',
                    dest='max_components',
                    help='Number of components to project and plot in 2D',
                    type=int,
                    default=2)

inputs.add_argument('-title',
                    action='store',
                    dest='title',
                    help='Title for graphs',
                    type=str,
                    required=True)

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output folder for text and png files',
                    type=str,
                    required=True)

# Parse into useful form
UserInput = parser.parse_args()

# Make output directory
out_dir=UserInput.out_dir
fname = './' + 'Pca_' + out_dir
if not os.path.exists(fname):
    os.mkdir(fname)

# Process trajectory
trajectory = TrajectoryReader.DCD(topology_path=UserInput.structure, trajectory_path=UserInput.trajectory).load()
trajectory = AtomSelection.Slice(trajectory=trajectory, atom_selection=UserInput.sel).select()
xyz = Featurizer.XYZ(trajectory=trajectory).extract()
del trajectory

# Reduce dimensions
projection, components, explained_variance = DimensionReduction.PCA(coordinates=xyz).reduce()

Saver.Array(array=projection, out_name=os.path.join('Pca_' + UserInput.out_dir, 'Projection.txt')).save()
Saver.Array(array=components, out_name=os.path.join('Pca_' + UserInput.out_dir , 'components.txt')).save()
Saver.Array(array=explained_variance, out_name=os.path.join('Pca_' + UserInput.out_dir, 'explained_variance.txt')).save()

# Plot graphs
for i in range(0, UserInput.max_components - 1):
    for j in range(i + 1, UserInput.max_components):
        energy, x_edges, y_edges = FreeEnergy.Rice(projection[:, i], projection[:, j]).calculate()
        Plotter.MeshContour(y=energy,
                            x_edges=x_edges,
                            y_edges=y_edges,
                            x_label='PC{0} ({1:.2f}%)'.format(i+1, explained_variance[i] * 100),
                            y_label='PC{0} ({1:.2f}%)'.format(j+1, explained_variance[j] * 100),
                            title=UserInput.title,
                            out_name=os.path.join('Pca_' + UserInput.out_dir, 'pc{0}_pc{1}.png'.format(i+1, j+1))
                            ).plot()

