
import argparse
import os
import numpy
from Analysis import AtomSelection, DimensionReduction, Featurizer, FreeEnergy, Plotter, Saver, TrajectoryReader

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
    description='Calculate, save and plot principal components and projections', add_help=False
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
inputs.add_argument('-sel',
                    action='store',
                    dest='sel',
                    help='Atom selection',
                    type=str,
                    default='name CA'
                    )
inputs.add_argument('-n',
                    action='store',
                    dest='max_components',
                    help='Number of components to project and plot in 2D',
                    type=int,
                    default=2
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

# Process trajectory
trajectory = TrajectoryReader.DCD(topology_path=UserInput.structure, trajectory_path=UserInput.trajectory).load()
trajectory = AtomSelection.Slice(trajectory=trajectory, atom_selection=UserInput.sel).select()
xyz = Featurizer.XYZ(trajectory=trajectory).extract()
del trajectory

projection, components, explained_variance = DimensionReduction.PCA(coordinates=xyz).reduce()

Saver.Array(array=projection, out_name=os.path.join(UserInput.out_name, 'projection.txt')).save()
Saver.Array(array=components, out_name=os.path.join(UserInput.out_name , 'components.txt')).save()
Saver.Array(array=explained_variance, out_name=os.path.join(UserInput.out_name, 'explained_variance.txt')).save()

for i in range(0, UserInput.max_components - 1):
    for j in range(i + 1, UserInput.max_components):
        energy, x_edges, y_edges = FreeEnergy.Rice(projection[:, i], projection[:, j]).calculate()
        Plotter.MeshContour(y=energy,
                           x_edges=x_edges,
                           y_edges=y_edges,
                           x_label='PCA_{0} ({1:.2f}%)'.format(i, explained_variance[i] * 100),
                           y_label='PCA_{0} ({1:.2f}%)'.format(j, explained_variance[j] * 100),
                           title='PCA_{0} and PCA_{1}'.format(i, j),
                           out_name=os.path.join(UserInput.out_name, 'PCA_{0}_PCA{1}.png'.format(i, j))
                           ).plot()
