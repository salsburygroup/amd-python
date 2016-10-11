#!/usr/env/ Python

# Ryan Melvin

import subprocess
import re
import math
import os
import hdbscan
import pandas
import numpy
import mdtraj
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial.distance import euclidean


# Initialize parser for user input
parser = argparse.ArgumentParser(description='Run intelligently sampled ACEMD using HDBSCAN', add_help=False)
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-d',
                    '--device',
                    action='store',
                    dest='device',
                    help='Which GPU card to use',
                    type=str,
                    required=True)
inputs.add_argument('-c',
                    '--config',
                    action='store',
                    dest='configuration_file',
                    help='Configuration file as if you were running a full sim',
                    type=str,
                    required=True)
inputs.add_argument('-sel',
                    action='store',
                    dest='sel',
                    help='Atom selection',
                    type=str,
                    default='not water')
UserInput = parser.parse_args()

with open(UserInput.configuration_file, "r") as input_conf:
    lines = input_conf.readlines()
with open(UserInput.configuration_file, "r") as input_conf:
    configuration_text = input_conf.read()
    
desired_steps = int(re.findall(r'run\s+(\d+)', configuration_text)[0])
topology = re.findall(r'structure\s+(.*)', configuration_text)[0]
pdb_file = re.findall(r'coordinates\s+(.*)', configuration_text)[0]
trajectory = re.findall(r'dcdfile\s+(.*)', configuration_text)[0]
sample_rate = int(re.findall(r'dcdfreq\s+(\d+)', configuration_text)[0])
desired_frames = int(math.ceil(desired_steps / sample_rate))
del configuration_text

with open(UserInput.configuration_file, "w") as input_conf:
    for line in lines:
        input_conf.write(re.sub(r'dcdfreq\s+\d+', r'dcdfreq\t1', line))
del lines

with open(UserInput.configuration_file, "r") as input_conf:
    lines = input_conf.readlines()

psf = mdtraj.load_topology(topology)
pdb = mdtraj.load(pdb_file)
atoms_of_interest = psf.select(UserInput.sel)
n_atoms_of_interest = len(atoms_of_interest)
reference_structure = pdb.atom_slice(atoms_of_interest)
not_water = psf.select("not water")

with mdtraj.formats.DCDTrajectoryFile('temp.dcd', 'w') as strided_trajectory:
    for output_frame in range(desired_frames):
        xyz = numpy.zeros((sample_rate, n_atoms_of_interest, 3), dtype='float64')
        if not os.path.exists(str(output_frame)):
            os.mkdir(str(output_frame))
        with open('input', "w+") as short_config:
            for line in lines:
                short_config.write(re.sub(r'run\s+\d+', r'run\t' + str(sample_rate * (output_frame + 1)), line))
        acemd_cmd = 'acemd --device ' + UserInput.device + ' input>&log_' + str(output_frame)
        process = subprocess.Popen(acemd_cmd, shell=True)
        process.wait()
        for frame in range(output_frame, output_frame + sample_rate):
            current_frame = mdtraj.load_frame(trajectory, frame, top=topology, atom_indices=atoms_of_interest)
            current_frame = current_frame.superpose(reference_structure)
            xyz[frame - output_frame, :, :] = current_frame.xyz
        xyz = xyz.astype('float64')
        clusterer = hdbscan.HDBSCAN(min_cluster_size=10)
        cluster_labels = clusterer.fit_predict(xyz)
        labeled_traj = pandas.DataFrame(columns=['frame', 'cluster'])
        labeled_traj['frame'] = numpy.arange(sample_rate)
        labeled_traj['cluster'] = cluster_labels
        for i in range(0, int(max(cluster_labels)) + 1):
            cluster_string = labeled_traj.loc[labeled_traj['cluster'] == i].frame.values
            cluster_coords = xyz[cluster_string]
            mean = cluster_coords.mean(axis=0)
            distance = [euclidean(row, mean) for row in cluster_coords]
            rep = cluster_string[numpy.argmin(distance)]
            mdtraj.load_frame(trajectory, rep, top=topology, atom_indices=not_water
                              )[rep].save_pdb(str(output_frame) + '/rep' + str(i) + '.pdb')
        del xyz
        plt.figure()
        plt.scatter(numpy.arange(sample_rate), cluster_labels, marker='+')
        plt.xlabel('Frame')
        plt.ylabel('Cluster')
        plt.title('HDBSCAN')
        plt.savefig(str(output_frame) + '/hdbscan_timeseries.png')
        plt.close()
        strided_trajectory.write(mdtraj.load_frame(trajectory, -1, top=topology).xyz)
        mdtraj.load_frame(trajectory, -1, top=topology).save_dcd(trajectory)
os.remove(trajectory)
os.rename('temp.dcd', trajectory)