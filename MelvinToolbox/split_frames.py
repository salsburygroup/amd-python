def split_DCD(structure, trajectory, output_prefix = 'Frame'):
    """
    ============
    Split frames
    ============

    Usage: %split(trajectory.dcd, output_prefix)
    Saves each frame in a DCD as a PDB with name output_prefix%s.pdb % frame#

    inputs:
        -structure: path to psf or pdb as string
        -trjectory: path to trajectory as string
        -output_prefix: prefix for output frames as string. If a prefix is not specified, 'Frame' will be used.

    :Author: Ryan Melvin
    :Date: Tue May 26 15:27:36 EDT 2015

    Sources: 
        - http://www.mdanalysis.org/MDAnalysisTutorial/writing.html

    Dependencies:
        - MDAnalysis

    """

    import MDAnalysis

    #Set up universe and selection
    u = MDAnalysis.Universe(structure, trajectory)
    sel = u.selectAtoms("all")
    
    #Loop over frames, saving each as a pdb
    for j in range(len(u.trajectory)):
        u.trajectory[j]
        name = output_prefix + str(j) + ".pdb" #prepare output name
        sel.write(name)



