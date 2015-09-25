#!/usr/bin/env python

#In its current implementation, this script could require memory up to three times the size of the input trajectory.

from __future__ import division
import mdtraj
import math
import os
import subprocess
import argparse
import numpy as np
from PIL import Image

# Find the helper file generate_shadow.vmd
dir = os.path.dirname(__file__)
cwd = os.getcwd()
shadow_helper = os.path.join(dir, 'generate_shadow.vmd')
median_helper = os.path.join(dir, 'generate_median.vmd')

parser = argparse.ArgumentParser(
        description = (
            'outputs a pdb with all frames within 1 sigma of the first frame in each cluster. Use the included vmd script to visualize these.' 
            ), 
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

        # Now, let's make some pretty pictures
        vmd_render_shadow_cmd = ('/Applications/VMD\ 1.9.2.app/Contents/vmd/vmd_MACOSXX86 ' 
            + directory +'/sigma.pdb -dispdev text -e ' 
            + shadow_helper + ' -args -first 1 -last ' + str(len(cluster)) + ' -outfile ' 
            + directory + '/shadow.dat')
        vmd_render_shadow=subprocess.call(vmd_render_shadow_cmd,shell=True)


        vmd_render_median_cmd = ('/Applications/VMD\ 1.9.2.app/Contents/vmd/vmd_MACOSXX86 ' 
                + directory +'/sigma.pdb -dispdev text -e ' 
                + median_helper + ' -args -outfile '
                +directory + '/median.dat'
                )
        vmd_render_median=subprocess.call(vmd_render_median_cmd,shell=True)
        
        tachyon_render_shadow_cmd = ('/Applications/VMD\ 1.9.2.app/Contents/vmd/tachyon_MACOSXX86'+
                ' -trans_vmd ' + directory + '/shadow.dat -o ' + directory + '/shadow.tga')
        tachyon_render_shadow=subprocess.call(tachyon_render_shadow_cmd,shell=True)

        tachyon_render_median_cmd = ('/Applications/VMD\ 1.9.2.app/Contents/vmd/tachyon_MACOSXX86'+
                ' -trans_vmd ' + directory + '/median.dat -o ' + directory + '/median.tga')
        tachyon_render_median=subprocess.call(tachyon_render_median_cmd,shell=True)

        #Let's get rid of the white pixels and convert the TGAs to PNGs
        median_img = Image.open(directory + '/median.tga')
        median_img = median_img.convert("RGBA")
        median_data = median_img.getdata()

        new_median_data = []
        for item in median_data:
            if item[0] == 89 and item[1] == 89 and item[2] == 89:
                new_median_data.append((89,89,89,0))
            else:
                new_median_data.append(item)

        median_img.putdata(new_median_data)
        median_img.save(directory + '/median.png',"PNG")

        shadow_img = Image.open(directory + '/shadow.tga')
        shadow_img = shadow_img.convert("RGBA")
        shadow_data = shadow_img.getdata()

        new_shadow_data = []
        for item in shadow_data:
            if item[0] == 255 and item[1] == 255 and item[2] == 255:
                new_shadow_data.append((255,255,255,0))
            else:
                new_shadow_data.append(item)

        shadow_img.putdata(new_shadow_data)
        shadow_img.save(directory + '/shadow.png',"PNG")

        #Now, let's layer them together
        layered_img = Image.alpha_composite(shadow_img,median_img)
        layered_img.save(directory + '/layered.png',"PNG")

        blended_img = Image.blend(median_img,shadow_img,0.5)
        blended_img.save(directory + '/blended.png',"PNG")

        # The user has specified a cutoff. Let's see if we've reached it.
        if cluster_number == UserInput.last:
            break
        else:
            cluster_number = cluster_number+1


