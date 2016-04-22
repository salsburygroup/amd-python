import mdtraj


class DCDReader:
    def __init__(self, dcd_path, topology_path):
        assert isinstance(dcd_path, str)
        self.dcd = dcd_path
        assert isinstance(topology_path, str)
        self.top = topology_path

    def load(self):
        trajectory = mdtraj.iterload(self.dcd, top=self.top)
        return trajectory