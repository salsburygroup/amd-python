#!/usr/bin/env python

#In its current implementation, this script could require memory up to three times the size of the input trajectory.

from __future__ import division
import mdtraj
import math
import os
import argparse
import numpy as np

parser = argparse.ArgumentParser(
        description = (
            'outputs a pdb with all frames within 1 sigma and another of those within 1 standard error of the first frame in each cluster. Use the included vmd script to visualize these.' 
            ), 
        add_help=False
        ) 
#List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-s', '-structure', action='store', dest='topology',help='Structure file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-t', '--traj', action='store', dest='dcd',help='Trajectory',type=str,required=True)
inputs.add_argument('-c', '--cluster', action='store', dest='cluster_data',help='cluster data file',type=str,required=True)


#Parse into useful form
UserInput=parser.parse_args()

trajectory = mdtraj.load(UserInput.dcd,top=UserInput.topology) #Load trajectory
cluster_number = 1 #We'll start counting clusters from 1
with open (UserInput.cluster_data) as file:
    for line in file:
        cluster = map(int,line.split( )) #From string to ints
        # Get frames in the user-specified distribution
        trajectory_subset = trajectory.slice(cluster)
        trajectory_subset.superpose(trajectory_subset,frame=0)

        # Now, we'll measure the RMSDs
        rmsd=mdtraj.rmsd(trajectory_subset,trajectory_subset,frame=0)

        # We'll save each cluster's statistics in its own folder
        directory = 'cluster' + str(cluster_number)
        os.makedirs(directory)

        # Save the frame used as the mean
        trajectory_subset[0].save(directory+'/mu.pdb')

        trajectory_subset = None #Won't be needing this again
        # Keep the RMSDs as a numpy array
        rmsd = np.array(rmsd)

        # Taking the cluster center to be the mean, calculate a modified standard deviation
        sum_squares = sum([x**2 for x in rmsd])
        sigma = math.sqrt(sum_squares/len(cluster))
        # And modified standard error
        SE = sigma/math.sqrt(len(cluster))

        # Now, let's use logical masks to get out those frames within sigma and SE
        sigma_mask = rmsd <= sigma
        SE_mask = rmsd <= SE
        sigma_frames = [frame for (frame,sigma_mask) in zip(cluster,sigma_mask) if sigma_mask]
        SE_frames = [frame for (frame,SE_mask) in zip(cluster,SE_mask) if SE_mask]
      

        with open(directory+'/sigma.txt','wb') as sigma_file:
            sigma_file.write(' '.join([str(i) for i in sigma_frames]))

        with open(directory+'/SE.txt','wb') as SE_file:
            SE_file.write(' '.join([str(i) for i in SE_frames]))

        # Let's also save the subsets of frames
        trajectory.slice(sigma_frames).save(directory+'/sigma.pdb')
        trajectory.slice(SE_frames).save(directory+'/SE.pdb')

        cluster_number = cluster_number+1
