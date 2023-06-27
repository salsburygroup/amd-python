import argparse
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

def calculate_rmsd(input_xtc, input_top):
    # Load the XTC file with the corresponding topology file
    trajectory = md.load_xtc(input_xtc, top=input_top)

    # Get the reference structure (the first frame)
    reference_frame = trajectory[0]

    # Calculate alpha carbon RMSD
    alpha_carbon_indices = trajectory.topology.select('name CA')
    alpha_carbon_rmsd = md.rmsd(trajectory, reference_frame, atom_indices=alpha_carbon_indices)
    alpha_carbon_rmsd_angstrom = alpha_carbon_rmsd * 10  # Convert to angstroms
    print("Alpha Carbon RMSD (in angstroms):", alpha_carbon_rmsd_angstrom)

    # Calculate all-atom RMSD
    all_atom_rmsd = md.rmsd(trajectory, reference_frame)
    all_atom_rmsd_angstrom = all_atom_rmsd * 10  # Convert to angstroms
    print("All-Atom RMSD (in angstroms):", all_atom_rmsd_angstrom)

    # Plot the RMSD values
    time = np.arange(trajectory.n_frames) * trajectory.timestep / 1000  # Time in picoseconds, converted to nanoseconds
    plt.plot(time, alpha_carbon_rmsd_angstrom, label='Alpha Carbon RMSD')
    plt.plot(time, all_atom_rmsd_angstrom, label='All-Atom RMSD')
    plt.xlabel('Time (ns)')
    plt.ylabel('RMSD (Å)')
    plt.legend()
    plt.title('RMSD of Trajectory')
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate the alpha carbon RMSD and all-atom RMSD for a given trajectory file based on the first structure.')
    parser.add_argument('input_xtc', type=str, help='Input XTC file')
    parser.add_argument('input_top', type=str, help='Input topology file (e.g., PDB)')

    args = parser.parse_args()

    calculate_rmsd(args.input_xtc, args.input_top)

