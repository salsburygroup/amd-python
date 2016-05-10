import matplotlib.pyplot
import numpy
import pyemma


class Plotter:
    def __init__(self, out_name):
        self.out_name = out_name

    def plot(self):
        raise NotImplementedError


class TimeSeries(Plotter):
    def __init__(self, out_name, labels):
        self.labels = labels
        super().__init__(out_name)

    def plot(self):
        frames = numpy.arange(self.labels.shape[0])
        matplotlib.pyplot.figure()
        matplotlib.pyplot.scatter(frames, self.labels, marker='+')
        matplotlib.pyplot.xlabel('Frame')
        matplotlib.pyplot.ylabel('Cluster')
        matplotlib.pyplot.title('Time Series')
        matplotlib.pyplot.savefig(self.out_name)
        matplotlib.pyplot.close()


class Scores(Plotter):
    def __init__(self, out_name, scores, parameter_list):
        self.scores = scores
        self.parameter_list = parameter_list
        super().__init__(out_name)

    def plot(self):
        matplotlib.pyplot.plot(self.parameter_list, self.scores)
        matplotlib.pyplot.xlabel('Parameter')
        matplotlib.pyplot.ylabel('CVI')
        matplotlib.pyplot.title('Scores')
        matplotlib.pyplot.savefig(self.out_name)
        matplotlib.pyplot.close()


class RateMatrix(Plotter):
    def __init__(self, out_name, msm):
        self.msm = msm
        super().__init__(out_name)

    def plot(self):
        matplotlib.pyplot.imshow(self.msm.transition_matrix,
                                 extent=(self.msm.active_set.min(), self.msm.active_set.max(),
                                         self.msm.active_set.min(), self.msm.active_set.max())
                                 )
        matplotlib.pyplot.colorbar()
        matplotlib.pyplot.xlabel('Transition to')
        matplotlib.pyplot.ylabel('Transition From')
        matplotlib.pyplot.title('Estimated Transition Matrix')
        matplotlib.pyplot.savefig(self.out_name)
        matplotlib.pyplot.clf()


class TransitionPath(Plotter):
    def __init__(self, out_name, msm):
        assert isinstance(msm, pyemma.msm.MaximumLikelihoodMSM)
        self.msm = msm
        super().__init__(out_name)

    def plot(self):
        pyemma.plots.plot_markov_model(self.msm)
        matplotlib.pyplot.title('Estimated Markov Chain')
        matplotlib.pyplot.savefig(self.out_name)
        matplotlib.pyplot.clf()
