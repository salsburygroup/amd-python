#from moleculekit.molecule import Molecule
from htmd.ui import *
mol = Molecule("ionized.pdb")
mol.read("output_stride100.dcd")
mol.wrap("protein") # or "protein and resid X to Y" if you want to wrap around specific protein residues
mol.write("output_wrapped_stride100.dcd")
