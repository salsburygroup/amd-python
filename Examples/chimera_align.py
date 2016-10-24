from MatchMaker import match, CP_SPECIFIC_SPECIFIC, GAP_OPEN, GAP_EXTEND, \
	defaults, MATRIX, ITER_CUTOFF

from chimera import openModels, Molecule

mol0 = chimera.openModels.open('5IIZ', type="PDB")
mol1 = chimera.openModels.open('5IM9', type="PDB")

s1, s2 = openModels.list(modelTypes=[Molecule])
c1, c2 = s1.sequence('A'), s2.sequence('A')

atoms1, atoms2, rmsd = match(CP_SPECIFIC_SPECIFIC, [(c1, c2)], defaults[MATRIX],
			"nw", defaults[GAP_OPEN], defaults[GAP_EXTEND],
			iterate=defaults[ITER_CUTOFF])[0]
print "RMSD:", rmsd
for a1, a2 in zip(atoms1, atoms2):
	print a1.residue, a2.residue, a1.xformCoord().distance(a2.xformCoord())
