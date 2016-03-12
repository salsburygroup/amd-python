import mdtraj


class DCDShaper:
    def __init__(self, trajectory, atom_selection):
        assert isinstance(trajectory, mdtraj.Trajectory)
        self.trajectory = trajectory
        assert isinstance(atom_selection, str)
        self.atom_selection = atom_selection

    def make_2d(self):
        sel = self.trajectory.topology.select(self.atom_selection)
        self.trajectory = self.trajectory.atom_slice(sel)
        frames = self.trajectory.n_frames
        atoms = self.trajectory.n_atoms
        trajectory_2d = self.trajectory.xyz.reshape((frames, atoms * 3))
        trajectory_2d = trajectory_2d.astype('float64')
        return trajectory_2d
