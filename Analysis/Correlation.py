import numpy
import mdtraj


class Correlation:
    def __init__(self, trajectory):
        assert isinstance(trajectory, mdtraj.Trajectory)
        self.trajectory = trajectory
        self.correlation_matrix = []

    def calculate(self):
        raise NotImplementedError


class Pearson(Correlation):
    def calculate(self):
        average = numpy.average(self.trajectory.xyz, axis=0)
        fluctuations = self.trajectory.xyz - average[numpy.newaxis, :]
        del average
        dots = numpy.zeros((self.trajectory.n_atoms, self.trajectory.n_atoms))
        for i in range(self.trajectory.n_frames):
            dot = numpy.dot(fluctuations[i, :, :], numpy.transpose(fluctuations[i, :, :]))
            dots = dots + dot
        del fluctuations
        dots = numpy.divide(dots, self.trajectory.n_frames)
        diagonal = numpy.diag(dots)
        normalization_matrix = numpy.outer(diagonal, diagonal)
        normalization_matrix = numpy.sqrt(normalization_matrix)
        self.correlation_matrix = numpy.divide(dots, normalization_matrix)
        return self.correlation_matrix


class MutualInformation(Correlation):
    def calculate(self):
        raise NotImplementedError
