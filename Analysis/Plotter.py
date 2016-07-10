import matplotlib
matplotlib.use('Agg') # For use on DEAC cluster
import matplotlib.pyplot
import numpy


class Plotter:
    def __init__(self, y, out_name, x_label=' ', y_label=' ', title=' '):
        self.y = y
        self.out_name = out_name
        self.x_label = x_label
        self.y_label = y_label
        self.title = title

    def plot(self):
        raise NotImplementedError


class Y(Plotter):
    def plot(self):
        matplotlib.pyplot.figure()
        matplotlib.pyplot.plot(self.y)
        matplotlib.pyplot.xlabel(self.x_label)
        matplotlib.pyplot.ylabel(self.y_label)
        matplotlib.pyplot.title(self.title)
        matplotlib.pyplot.savefig(self.out_name)
        matplotlib.pyplot.clf()


class XY(Plotter):
    def __init__(self, x, y, out_name, x_label=' ', y_label=' ', title=' '):
        self.x = x
        super().__init__(y, out_name, x_label, y_label, title)

    def plot(self):
        matplotlib.pyplot.figure()
        matplotlib.pyplot.plot(self.x, self.y)
        matplotlib.pyplot.xlabel(self.x_label)
        matplotlib.pyplot.ylabel(self.y_label)
        matplotlib.pyplot.title(self.title)
        matplotlib.pyplot.savefig(self.out_name)
        matplotlib.pyplot.clf()


class SimplePColor(Plotter):
    def plot(self):
        matplotlib.pyplot.pcolor(self.y, cmap='jet')
        matplotlib.pyplot.colorbar()
        matplotlib.pyplot.xlabel(self.x_label)
        matplotlib.pyplot.ylabel(self.y_label)
        matplotlib.pyplot.title(self.title)
        matplotlib.pyplot.savefig(self.out_name)
        matplotlib.pyplot.clf()


class MeshPColor(Plotter):
    def __init__(self, y, x_edges, y_edges, out_name, x_label=' ', y_label=' ', title=' '):
        self.x_edges = x_edges
        self.y_edges = y_edges
        super().__init__(y, out_name, x_label, y_label, title)

    def plot(self):
        masked_y = numpy.ma.masked_where(numpy.isnan(self.y), self.y)
        matplotlib.pyplot.figure()
        matplotlib.pyplot.pcolormesh(self.x_edges, self.y_edges, masked_y.T)
        matplotlib.pyplot.colorbar()
        matplotlib.pyplot.xlabel(self.x_label)
        matplotlib.pyplot.ylabel(self.y_label)
        matplotlib.pyplot.title(self.title)
        matplotlib.pyplot.savefig(self.out_name)
        matplotlib.pyplot.clf()
