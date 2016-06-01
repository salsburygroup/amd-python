# Visual Statistics
The two python scripts here generate an image of a biopolymer with uncertainty in position shown as shadow.

The script vsigma.py takes a structure file, a trajectory and clustering data where a row contains the frames in a cluster, with the first entry in the row being the median structure in the cluster.

The script vdistribution.py takes two input PDBs -- the median structure as one PDB and the part of the distribution to be displayed as shadow as the second PDB.
