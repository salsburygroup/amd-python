import matplotlib.pyplot
import numpy


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
        matplotlib.pyplot.clf()


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
        matplotlib.pyplot.clf()