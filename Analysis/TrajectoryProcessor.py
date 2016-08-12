import mdtraj


class Processor:
    def __init__(self, trajectory, atom_selection):
        assert isinstance(trajectory, mdtraj.Trajectory)
        self.trajectory = trajectory
        self.sel = atom_selection

    def process(self):
        raise NotImplementedError


class AtomIndexer(Processor):
    # From Oliver Schillinger's mdtraj tools https://github.com/schilli/Tools
    def process(self):
        self.atom_indices = self.trajectory.topology.select(self.sel)
        return self.atom_indices


class Stripper(Processor):
    def process(self):
        raise NotImplementedError


class Strider(Processor):
    def process(self):
        raise NotImplementedError


class Aligner(Processor):
    def process(self):
        indices = AtomIndexer(self.trajectory, self.sel).process()
        aligned_trajectory = self.trajectory.superpose(self.trajectory, atom_indices=indices)
        return aligned_trajectory


class Wrapper(Processor):
    def process(self):
        raise NotImplementedError


class Unwrapper(Processor):
    def process(self):
        raise NotImplementedError


