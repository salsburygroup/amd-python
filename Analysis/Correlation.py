import numpy
import mdtraj
import hdbscan
import os
import pandas
import signal
import subprocess


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


class Propagator(Correlation):
    def __init__(self, trajectory, tau):
        self.tau = tau
        super().__init__(trajectory)

    def calculate(self):
        delta_sum = numpy.zeros([self.trajectory.topology.n_atoms, 3], dtype=float)
        dot_sum = numpy.zeros([self.trajectory.topology.n_atoms, self.trajectory.topology.n_atoms], dtype=float)
        for frame in numpy.arange(self.trajectory.n_frames - self.tau):
            delta_temp = self.trajectory.xyz[frame] - self.trajectory.xyz[frame + self.tau]
            delta_sum = delta_sum + delta_temp
            dot_sum = dot_sum + numpy.inner(delta_temp, delta_temp)
        average_delta = delta_sum / (self.trajectory.n_frames - self.tau)
        average_dot = dot_sum / (self.trajectory.n_frames - self.tau)
        dot_average_delta = numpy.inner(average_delta, average_delta)
        return average_dot, average_delta, dot_average_delta


class Clustering(Correlation):
    def __init__(self, trajectory, input_type='correlation', minimum_membership=None, correlation_matrix=None):
        self.minimum_membership = minimum_membership
        input_options = ('correlation', 'similarity', 'dots')
        if input_type not in input_options:
            raise ValueError('input types are correlation, similarity or dots')
        self.input_type = input_type
        self.correlation_matrix = correlation_matrix
        self.labels = []
        super().__init__(trajectory)

    def calculate(self):
        try:
            self.correlation_matrix[0]
        except IndexError:
            self.correlation_matrix = Pearson(self.trajectory).calculate()

        number_residues = self.trajectory.topology.n_atoms
        three_percent = int(numpy.ceil(number_residues * 0.03))
        if self.minimum_membership:
            min_cluster_size = self.minimum_membership
        elif three_percent >= 2:
            min_cluster_size = three_percent
        else:
            min_cluster_size = 2

        if self.input_type == 'similarity':
            clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, metric='precomputed')
            distance = 1 - numpy.abs(self.correlation_matrix)
            labels = clusterer.fit_predict(distance)
        elif self.input_type == 'dots':
            clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, metric='precomputed')
            dots = self.correlation_matrix.dot(self.correlation_matrix.T)
            distance = 1 - dots / numpy.max(numpy.abs(dots))
            labels = clusterer.fit_predict(distance)
        else:
            clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size)
            labels = clusterer.fit_predict(self.correlation_matrix)

        return labels

    @staticmethod
    def visualize(labels, pdb_file, out_name):
        num_residues = len(labels)
        df = pandas.DataFrame(columns=['residue', 'cluster'])
        df['residue'] = numpy.arange(num_residues)
        df['cluster'] = labels
        clusters = numpy.unique(df.loc[df['cluster'] != -1].cluster.values)

        with open(out_name + '_image_script.vmd', 'w+') as vmd_file:
            vmd_file.write(
                'mol new ' + pdb_file + '\n' + 'mol delrep 0 top\n' + 'mol representation NewCartoon\n'
                + 'mol selection {all}\n' + 'mol material Ghost\n' + 'mol addrep top\n' + 'mol representation Bonds\n'
                + 'mol material AOChalky\n' + ' display ambientocclusion on\n')
            for cluster in clusters:
                cluster_string = ' '.join(
                    ['%d' % num for num in df.loc[df['cluster'] == cluster].residue.values]
                )
                vmd_file.write(
                    'mol color ColorID ' + str(cluster) + '\n' + 'mol selection {residue ' + cluster_string + '}\n'
                    + 'mol addrep top\n'
                )
            vmd_file.write(
                'display resize 1920 1080\n' + 'display resetview\n' + 'render TachyonInternal ' + out_name
                + '_render.tga\n' + 'exit'
            )

        # Posix-compliant way to exit shell
        signal.signal(signal.SIGTTOU, signal.SIG_IGN)
        # Now, let's make some pretty pictures
        vmd_render_cmd = (
            'vmd '
            + ' -dispdev text -e '
            + out_name + '_image_script.vmd'
        )
        subprocess.call([os.getenv('SHELL'), '-i', '-c', vmd_render_cmd])
        os.tcsetpgrp(0, os.getpgrp())
