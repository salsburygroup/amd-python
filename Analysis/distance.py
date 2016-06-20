import numpy
import mdtraj
from .AtomSelection import Slice


class Distance:
    def __init__(self, trajectory, atom_selection):
        assert isinstance(trajectory, mdtraj.Trajectory)
        self.trajectory = trajectory
        self.sel = atom_selection

    def calculate(self):
        raise NotImplementedError


class RMSF(Distance):
    def calculate(self):
        sub_trajectory = Slice(trajectory=self.trajectory, atom_selection=self.sel).select()
        assert isinstance(sub_trajectory, mdtraj.Trajectory)
        reference_positions = sub_trajectory.xyz.mean(0)
        difference = sub_trajectory.xyz - reference_positions
        sum_squares = numpy.sum(numpy.sum(numpy.square(difference), axis=2), axis=0)
        # Code at https://github.com/schilli/Tools/blob/master/rmsf.py uses 3*n_frames. Seems wrong...
        rmsf = (sum_squares/sub_trajectory.n_frames)**0.5
        return rmsf
