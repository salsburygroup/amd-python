#!/usr/bin/env python

#In its current implementation, this script could require memory up to three times the size of the input trajectory.

from __future__ import division
import mdtraj
import math
import os
import argparse
import numpy as np
from PIL import Image

# Find the helper file generate_shadow.vmd
dir = os.path.dirname(__file__)
cwd = os.getcwd()
shadow_helper = os.path.join(dir, 'generate_shadow.vmd')
median_helper = os.path.join(dir, 'generate_middle.vmd')
tachyon = os.path.join(dir, 'tachyon')

parser = argparse.ArgumentParser(
        description = (
            'outputs a pdb with all frames within 1 sigma of the first frame in each cluster.'), 
        add_help=False
        ) 
#List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-s', '--structure', action='store', dest='topology',help='Structure file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-t', '--traj', action='store', dest='dcd',help='Trajectory',type=str,required=True)
inputs.add_argument('-c', '--cluster', action='store', dest='cluster_data',help='cluster data file',type=str,required=True)
inputs.add_argument('-l', '--last', action='store', dest='last',help='last cluster for which to generate files',type=int,required=True)

#Parse into useful form
UserInput=parser.parse_args()

trajectory = mdtraj.load(UserInput.dcd,top=UserInput.topology) #Load trajectory
cluster_number = 1 #We'll start counting clusters from 1
with open (UserInput.cluster_data) as file:
    for line in file:
        cluster = map(int,line.split( )) #From string to ints
        # Get framesusin the user-specified distribution
        trajectory_subset = trajectory.slice(cluster)
        trajectory_subset.superpose(trajectory_subset,frame=0) #Subtracting off the median

        # Now, we'll measure the RMSDs
        rmsd=mdtraj.rmsd(trajectory_subset,trajectory_subset,frame=0)

        # We'll save each cluster's statistics in its own folder
        directory = cwd + '/cluster' + str(cluster_number)
        os.makedirs(directory)

        # Save the frame used as the mean
        trajectory_subset[0].save(directory+'/mu.pdb')

        #trajectory_subset = None #Won't be needing this again
        # Keep the RMSDs as a numpy array
        rmsd = np.array(rmsd)

        # Taking the cluster median to be the mean, calculate a modified standard deviation
        sum_squares = sum([x**2 for x in rmsd])
        sigma = math.sqrt(sum_squares/len(cluster))

        # Now, let's use logical masks to get out those frames within sigma 
        sigma_mask_txt = rmsd <= sigma
        sigma_mask = [i for i in range(len(rmsd)) if rmsd[i] <= sigma]
        sigma_frames_txt = [frame for (frame,sigma_mask_txt) in zip(cluster,sigma_mask_txt) if sigma_mask_txt]
      

        with open(directory+'/sigma.txt','wb') as sigma_file:
            sigma_file.write(' '.join([str(i) for i in sigma_frames_txt]))

        # Let's also save the subsets of frames
        trajectory_subset.slice(sigma_mask).save(directory+'/sigma.pdb')
        
