import mdtraj


class Reader:
    def __init__(self, trajectory_path):
        self.trajectory = trajectory_path

    def load(self):
        raise NotImplementedError


class DCD(Reader):
    def __init__(self, trajectory_path, topology_path):
        self.topology = topology_path
        super().__init__(trajectory_path)

    def load(self):
        trajectory = mdtraj.load(self.trajectory, top=self.topology)
        return trajectory


class BigDCD(Reader):
    def load(self):
        raise NotImplementedError