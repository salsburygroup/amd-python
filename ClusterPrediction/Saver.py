import numpy


class Saver:
    def __init__(self, out_name):
        self.out_name = out_name

    def save(self):
        raise NotImplementedError


class TimeSeries(Saver):
    def __init__(self, out_name, labels):
        self.labels = labels
        super().__init__(out_name)

    def save(self):
        numpy.savetxt(self.out_name, self.labels)


class Scores(Saver):
    def __init__(self, out_name, scores):
        self.scores = scores
        super().__init__(out_name)

    def save(self):
        numpy.savetxt(self.out_name, self.scores)