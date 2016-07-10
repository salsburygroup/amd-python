#!/usr/bin/env python

# Example call
# python /Users/melvrl13/Documents/AMD/AMD-PYTHON/bin/vmd_qt_representatives.py -s /Volumes/RyanMdata-1/FUMP10/Folding/weightedSims3200/hairpin.pdb -t /Volumes/RyanMdata-1/FUMP10/Folding/weightedSims3200/weighted3200.dcd -c /Volumes/RyanMdata-1/FUMP10/Folding/weightedSims3200/clustering/cluster5.dat -o /Volumes/RyanMdata-1/FUMP10/Folding/weightedSims3200/repTest

import mdtraj
import os
import argparse

cwd = os.getcwd()

parser = argparse.ArgumentParser(
        description = (
            'Extracts cluster representatives as PDBs given input trajectory and clustering data with rows as clusters.'
            ), 
        add_help=False
        ) 
#List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-s', '--structure', action='store', dest='topology',help='Structure file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-t', '--traj', action='store', dest='dcd',help='Trajectory',type=str,required=True)
inputs.add_argument('-c', '--cluster', action='store', dest='cluster_data',help='cluster data file',type=str,required=True)
inputs.add_argument('-o', '--outdir', action='store', dest='cwd',help='output directory',type=str, default=cwd)

#Parse into useful form
UserInput=parser.parse_args()

trajectory = mdtraj.load(UserInput.dcd,top=UserInput.topology) #Load trajectory
cluster_number = 1 #We'll start counting clusters from 1
cwd = UserInput.cwd

with open (UserInput.cluster_data) as file:
    for line in file:
        cluster = list(map(int,line.split( ))) #From string to ints
        trajectory_subset = trajectory.slice(cluster)
        trajectory_subset.superpose(trajectory_subset,frame=0) #Subtracting off the median
        trajectory_subset[0].save(UserInput.cwd + '/rep{0}.pdb'.format(cluster_number))
        cluster_number += 1