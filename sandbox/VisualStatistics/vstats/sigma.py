def sigma(cluster, trajectory, cluster_number, directory):
        # Get frames in the user-specified distribution
        trajectory_subset = trajectory.slice(cluster)
        trajectory_subset.superpose(trajectory_subset,frame=0) #Subtracting off the median

        # Now, we'll measure the RMSDs
        rmsd=mdtraj.rmsd(trajectory_subset,trajectory_subset,frame=0)

        # We'll save each cluster's statistics in its own folder
        directory = cwd + '/cluster' + str(cluster_number)
        os.makedirs(directory)

        # Save the frame used as the mean
        trajectory_subset[0].save(directory+'/mu.pdb')

        trajectory_subset = None #Won't be needing this again
        # Keep the RMSDs as a numpy array
        rmsd = np.array(rmsd)

        # Taking the cluster median to be the mean, calculate a modified standard deviation
        sum_squares = sum([x**2 for x in rmsd])
        sigma = math.sqrt(sum_squares/len(cluster))

        # Now, let's use logical masks to get out those frames within sigma 
        sigma_mask = rmsd <= sigma
        sigma_frames = [frame for (frame,sigma_mask) in zip(cluster,sigma_mask) if sigma_mask]
      
        with open(directory+'/sigma.txt','wb') as sigma_file:
            sigma_file.write(' '.join([str(i) for i in sigma_frames]))

        # Let's also save the subsets of frames
        trajectory.slice(sigma_frames).save(directory+'/sigma.pdb')


