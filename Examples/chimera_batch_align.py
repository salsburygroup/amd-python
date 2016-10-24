import chimera

def findMAV():
	# locate and return the (presumably only) running instance
	# of MultAlign Viewer
	from MultAlignViewer.MAViewer import MAViewer
	from chimera.extension import manager
	mavs = [inst for inst in manager.instances
					if isinstance(inst, MAViewer)]
	if len(mavs) != 1:
		raise AssertionError("not exactly one MAV")
	return mavs[0]


# read a file named "structList" whose lines have 3 fields:  two PDB files
# and an output file name for the alignment
for line in open("structList", "r"):
	pdb1, pdb2, output = line.strip().split()

	# open the pdb files
	mol1 = chimera.openModels.open(pdb1, type="PDB")[0]
	mol2 = chimera.openModels.open(pdb2, type="PDB")[0]

	# use MatchMaker to get them superimposed
	from MatchMaker import match, CP_BEST, defaultMatrix, \
				defaultAlignAlgorithm, defaultGapOpen, \
				defaultGapExtend, defaultIterateCutoff
	# match the best-aligning chain of pdb1 onto pdb2, specifying
	# the similarity matrix, alignment algorithm, and the gap
	# open/extend penalties
	atoms1, atoms2, rmsd = match(CP_BEST, (mol1, [mol2]), defaultMatrix,
		defaultAlignAlgorithm, defaultGapOpen, defaultGapExtend,
		iterate=defaultIterateCutoff, showAlignment=True)[0]

	# find the Multalign Viewer that is now showing
	mav = findMAV()

	# write out the alignment in MSF format
	# other format choices:  Clustal ALN, Aligned FASTA, Aligned NBRF/PIR,
	#	Pfam, GCG RSF, and Stockholm
	from MultAlignViewer.output import saverFunc
	out = open(output, "w")
	saverFunc["MSF"](out, mav, mav.seqs, mav.fileMarkups)
	out.close()

	# write out the RMSD
	out = open(output + ".rmsd", "w")
	print>>out, "%.3f" % rmsd
	out.close()

	# again quit the MAV
	mav.Quit()

	# close the PDBs
	chimera.openModels.close([mol1, mol2])
raise chimera.ChimeraSystemExit("script finished\n")
