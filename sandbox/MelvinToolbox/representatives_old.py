def pdb(cluster, structure, trajectory, output_prefix = 'Frame'):

    """
    ===================
    Representatives.pdb
    ===================

    Usage: %representatives.pdb(structure, trajectory.dcd, frame_list, output_prefix)
    Saves each frame in list from DCD as a PDB with name output_prefix%s.pdb % frame#

    inputs:
        -cluster: space-separated text file of cluster data
        -structure: path to psf or pdb as string
        -trjectory: path to trajectory as string
        -frame_list: list of desired frames
        -output_prefix: prefix for output frames as string. If a prefix is not specified, 'Frame' will be used.

    :Author: Ryan Melvin
    :Date: Mon Jun  1 14:06:10 EDT 2015
    Sources: 

    Dependencies:
        - MDAnalysis
        - split_frames

    """
    from split_frames import split_list

    #Initialize empty list of cluster representatives
    representatives = []
    
    #Get list of representative frame number (should be in 0-based indexing)
    with open(cluster) as f:
        for line in f:
            representatives.append(line.split(None,1)[0])

    #Change type str to type int
    representatives = [int(x) for x in representatives]

    #Pass list and trajectory information to split_list
    split_list(structure, trajectory, representatives, output_prefix)
            
