def two_atoms(structure, trajectory, atom1, atom2, outname):
    """
    =========
    Two Atoms
    =========

    Usage: %two_atoms(structure.psf trajectory.dcd, output_name.txt)
    Time series of distance between two atoms

    inputs:
        -structure: path to psf or pdb as string
        -trjectory: path to trajectory as string
        -atom1: 1-based index number of atom as int
        -atom2: 1-based index number of atom as int
        -output_prefix: name for output file as string 

    :Author: Ryan Melvin
    :Date: Thu May 28 15:59:16 EDT 2015

    Sources: 
        - http://pythonhosted.org/MDAnalysis/documentation_pages/analysis/distances.html

    Dependencies:
        - MDAnalysis
        - Pandas
    """
    
    from MDAnalysis import Universe
    from MDAnalysis.analysis.distances import dist
    from pandas import DataFrame

    #Initialize variables
    d = [] #Temporary storage space for distance vector of the form (residue id of atom1, residue id of atom2, distance between atom1 and atom2)
    ts = DataFrame(columns=['R']) #Where the distance for each frame will be stored

    #Set up the universe
    u = Universe(structure,trajectory)

    #Select the two atoms
    ag1 = u.selectAtoms("bynum %s" % atom1)
    ag2 = u.selectAtoms("bynum %s" % atom2)

    #Find the distance for each frame
    for i in range(len(u.trajectory)):
        u.trajectory[i] #Go to current frame
        d = dist(ag1,ag2)
        ts.loc[i]=d[2] #Store only the distance

    #Save time series as csv in text file
    ts.to_csv(outname, index=False)
