import matplotlib.pyplot


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
        super().__init__(y, x_label, y_label, title, out_name)

    def plot(self):
        matplotlib.pyplot.figure()
        matplotlib.pyplot.plot(self.x, self.y)
        matplotlib.pyplot.xlabel(self.x_label)
        matplotlib.pyplot.ylabel(self.y_label)
        matplotlib.pyplot.title(self.title)
        matplotlib.pyplot.savefig(self.out_name)
        matplotlib.pyplot.clf()


class HeatMap(Plotter):
    def plot(self):
        matplotlib.pyplot.pcolor(self.y, cmap='jet')
        matplotlib.pyplot.colorbar()
        matplotlib.pyplot.xlabel(self.x_label)
        matplotlib.pyplot.ylabel(self.y_label)
        matplotlib.pyplot.title(self.title)
        matplotlib.pyplot.savefig(self.out_name)