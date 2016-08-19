import mdtraj
import numpy
from tempfile import mkdtemp
import os.path as path


class CorrelationPropagator:
    def __init__(self, dcd_path, top_path, atom_selection, tau=1):
        self.dcd_path = dcd_path
        self.top_path = top_path
        assert isinstance(atom_selection, str)
        self.atom_selection = atom_selection
        self.tau = tau
        self.average_delta = []
        self.average_dot = []
        self.dot_average_delta = []

    def matrix(self):
        trajectory = mdtraj.load(self.dcd_path, top=self.top_path)
        deltas = numpy.empty([trajectory.n_frames-self.tau, trajectory.topology.n_atoms, 3], dtype=float)
        temp_file = path.join(mkdtemp(), 'newfile.dat')
        # dots = numpy.empty([trajectory.n_frames-self.tau, trajectory.topology.n_atoms, trajectory.topology.n_atoms],
        #                   dtype=float
        #                   )
        dots = numpy.memmap(temp_file, dtype='float', mode='w+',
                            shape=(trajectory.n_frames-self.tau,
                                   trajectory.topology.n_atoms,
                                   trajectory.topology.n_atoms
                                   )
                            )

        for frame in numpy.arange(trajectory.n_frames - self.tau):
            deltas[frame, :, :] = trajectory.xyz[frame] - trajectory.xyz[frame+self.tau]
            dots[frame, :, :] = numpy.inner(deltas[frame, :, :], deltas[frame, :, :])
        self.average_delta = deltas.mean(axis=0)
        self.average_dot = dots.mean(axis=0)
        self.dot_average_delta = numpy.inner(self.average_delta, self.average_delta)
        return self.average_dot


if __name__ == "__main__":
    import argparse
    import os
    import matplotlib.pyplot as plt
    # Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
    parser = argparse.ArgumentParser(description='Calculate Propagating Correlation Matrix', add_help=False)

    # List all possible user input
    inputs = parser.add_argument_group('Input arguments')
    inputs.add_argument('-h', '--help', action='help')
    inputs.add_argument('-top', action='store', dest='structure',
                        help='Structure file corresponding to trajectory', type=str, required=True)
    inputs.add_argument('-traj', action='store', dest='trajectory', help='Trajectory', type=str, required=True)
    inputs.add_argument('-sel', action='store', dest='sel', help='Atom selection', type=str, default='not element H')
    inputs.add_argument('-tau', action='store', dest='tau', help='lag time', type=int, default=1)
    inputs.add_argument('-o', action='store', dest='out_name', help='Output directory', type=str, required=True)

    # Parse into useful form
    UserInput = parser.parse_args()

    # Execute calculation
    cp = CorrelationPropagator(UserInput.trajectory, UserInput.structure, UserInput.sel,UserInput.tau)
    average_dot = cp.matrix()

    # Save text results
    numpy.savetxt(os.path.join(UserInput.out_name, 'average_dot.txt'), average_dot)
    numpy.savetxt(os.path.join(UserInput.out_name, 'average_delta.txt'), cp.average_delta)
    numpy.savetxt(os.path.join(UserInput.out_name, 'dot_average_delta.txt'), cp.dot_average_delta)

    # Save pretty pictures
    plt.matshow(average_dot)
    plt.xlabel('Atom')
    plt.ylabel('Atom')
    plt.title('Average dot product')
    plt.colorbar()
    plt.savefig(os.path.join(UserInput.out_name, 'average_dot.png'))
    plt.close()

    plt.figure()
    plt.matshow(cp.dot_average_delta)
    plt.xlabel('Atom')
    plt.ylabel('Atom')
    plt.title('Dot product of average delta')
    plt.colorbar(format='%.0e')
    plt.savefig(os.path.join(UserInput.out_name, 'dot_average_delta.png'))
    plt.close()
