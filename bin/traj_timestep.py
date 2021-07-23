import mdtraj as md

t=md.load('output.dcd',top='ionized.pdb')
print(t.timestep)
