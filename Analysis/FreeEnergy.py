import numpy
from . import Histogrammer

class TwoDimensions:
    def __init__(self, vector1, vector2):
        self.vector1 = vector1
        self.vector2 = vector2

    def calculate(self):
        raise NotImplementedError


class FreedmanDiaconis(TwoDimensions):
    def calculate(self):
        histogram, x_edges, y_edges = Histogrammer.FreedmanDiaconis(
            vector1=self.vector1,
            vector2=self.vector2
        ).bin()
        histogram[histogram == 0] = numpy.nan
        energy_landscape = -0.6 * numpy.log(histogram)
        return energy_landscape, x_edges, y_edges
