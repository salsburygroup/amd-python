import mdtraj
from .TrajectoryReader import DCDReader


class Processor:
    def __init__(self, dcd_path, top_path, atom_selection):
        self.dcd = dcd_path
        self.top = top_path
        self.sel = atom_selection

    def process(self):
        raise NotImplementedError


class AtomIndexer(Processor):
    # From Oliver Schillinger's mdtraj tools https://github.com/schilli/Tools
    def process(self):
        trajectory = DCDReader(dcd_path=self.dcd, topology_path=self.top)
        assert isinstance(trajectory, mdtraj.Trajectory)
        atom_indices = trajectory.topology.select(self.sel)
        return atom_indices


class Stripper(Processor):
    def process(self):
        raise NotImplementedError


class Strider(Processor):
    def process(self):
        raise NotImplementedError


class Aligner(Processor):
    def process(self):
        raise NotImplementedError
