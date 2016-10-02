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
with open(UserInput.configuration_file, "w") as input_conf:
    for line in lines:
        input_conf.write(re.sub(r'dcdfreq\s+\d+', r'dcdfreq\t1', line))
del lines

desired_steps = int(re.findall(r'run\s+(\d+)', configuration_text)[0])
topology = re.findall(r'structure\s+(.*)', configuration_text)[0]
trajectory = re.findall(r'dcdfile\s+(.*)', configuration_text)[0]
tenth_nanoseconds = int(math.ceil(desired_steps / 25000))
del configuration_text

with open(UserInput.configuration_file, "r") as input_conf:
    lines = input_conf.readlines()

for p1_ns in range(tenth_nanoseconds):
    os.mkdir(str(p1_ns))
    with open('input', "w+") as short_config:
        for line in lines:
            short_config.write(re.sub(r'run\s+\d+', r'run\t' + str(25000 * (p1_ns + 1)), line))
    acemd_cmd = 'acemd --device ' + UserInput.device + ' input>&log_' + str(p1_ns)
    process = subprocess.Popen(acemd_cmd, shell=True)
    process.wait()
    t = mdtraj.load(trajectory, top=topology)
    sel = t.topology.select(UserInput.sel)
    t_slice = t.atom_slice(sel)
    t_slice = t_slice.superpose(t_slice)
    t_slice = t_slice[p1_ns:25000 + p1_ns - 1]
    temp = t_slice.xyz
    frames = t_slice.xyz.shape[0]
    atoms = t_slice.xyz.shape[1]
    data = temp.reshape((frames, atoms * 3))
    data = data.astype('float64')
    del temp
    clusterer = hdbscan.HDBSCAN(min_cluster_size=10)
    cluster_labels = clusterer.fit_predict(data)
    num_frames = len(cluster_labels)
    labeled_traj = pandas.DataFrame(columns=['frame', 'cluster'])
    labeled_traj['frame'] = numpy.arange(num_frames)
    labeled_traj['cluster'] = cluster_labels
    for i in range(0, int(max(cluster_labels)) + 1):
        cluster_string = labeled_traj.loc[labeled_traj['cluster'] == i].frame.values
        cluster_coords = data[cluster_string]
        mean = cluster_coords.mean(axis=0)
        distance = [euclidean(row, mean) for row in cluster_coords]
        rep = cluster_string[numpy.argmin(distance)]
        t_slice[rep].save_pdb(str(p1_ns) + '/rep' + str(i) + '.pdb')
    plt.figure()
    plt.scatter(numpy.arange(frames), cluster_labels, marker='+')
    plt.xlabel('Frame')
    plt.ylabel('Cluster')
    plt.title('HDBSCAN')
    plt.savefig(str(p1_ns) + '/hdbscan_timeseries.png')
    plt.close()
    del t_slice
    t[-1:-1 * (p1_ns + 1)].save_dcd(trajectory)
    del t
