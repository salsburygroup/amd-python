#####
#pca2.py   
#
#Ryan Melvin
#####
#Credit:
#If you use this script, cite
#TODO:
#####



#Outputs: aligned dcd

#Example Call
#python ~/Documents/AMD/AMD-PYTHON/pca.py '/Volumes/RyanMdata/F10/Folding/weightedSims3200/f10.psf' '/Volumes/RyanMdata/F10/Folding/weightedSims3200/weighted3200.dcd' '/Users/melvrl13/Desktop/testPCA/test/'

def xyz(topology,trajectory,outname):
    import mdtraj as md
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    from numpy import savetxt

    # Load and align trajectory
    traj = md.load(trajectory,top=topology)
    traj.superpose(traj,0) #Aligns trajectory

    # Create a 2-conponent model
    pca1=PCA(n_components=2) #initialize in memory
    reduced_cartesian = pca1.fit_transform(traj.xyz.reshape(traj.n_frames, traj.n_atoms * 3))

    # Prepare a plot
    plt.figure()
    plt.scatter(reduced_cartesian[:, 0], reduced_cartesian[:,1], marker='x', c=traj.time)
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('Cartesian coordinate PCA: ' + trajectory)
    cbar = plt.colorbar()
    cbar.set_label('Time')

    plt.savefig(outname+'.png') #save plot
    savetxt(outname+'.csv',reduced_cartesian,delimiter=',') #save projected PCA data

def pairwise(topology,trajectory,outname):
    import mdtraj as md
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    from itertools import combinations
    from numpy import savetxt

    # Load and align trajectory
    traj = md.load(trajectory,top=topology)

    # Get pairwise atom distances
    atom_pairs=list(combinations(range(traj.n_atoms),2))
    pairwise_distances= md.geometry.compute_distances(traj,atom_pairs)
    
    # Perform PCA on pairwise distances
    pca2 = PCA(n_components=2)
    reduced_distances = pca2.fit_transform(pairwise_distances)

    # Prepare a plot
    plt.figure()
    plt.scatter(reduced_distances[:, 0], reduced_distances[:,1], marker='x', c=traj.time)
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('Pairwise distance PCA: ' + trajectory)
    cbar = plt.colorbar()
    cbar.set_label('Time [ps]')

    plt.savefig(outname+'.png') #save plot
    savetxt(outname+'.csv',reduced_distances,delimiter=',') #save projected PCA data

if __name__ == "__main__":
    import sys
    xyz(sys.argv[1],sys.argv[2],sys.argv[3]+'_xyz')
    pairwise(sys.argv[1],sys.argv[2],sys.argv[3]+'_pairwise')
