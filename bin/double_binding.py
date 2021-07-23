#! /usr/bin/env python
from scipy.spatial import distance
import numpy as np
import MDAnalysis
import argparse
import os

# Find nth smallest number in a list 
def find_sub_min(brr, n):
    arr = brr.copy()
    for i in range(n-1):
        arr_ = arr
        arr_[np.argmin(arr_)] = np.max(arr)
        arr = arr_
    #print("# arr中最小的数为{}，位于第{}位".format(np.min(arr_), np.argmin(arr_)+1))
    return np.argmin(arr_)

# Calculate euclidean distance between two points
def euclidean_distance(a, b):
    d=((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)**(1/2)
    return d

print('Make sure to use the absolute path for output directory(-o)')
# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Compute the vector difference between Na+ binding loop and Na+', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-s',
                    action='store',
                    dest='structure',
                    help='Structure file corresponding to trajectory',
                    type=str,
                    required=True)

inputs.add_argument('-t',
                    action='store',
                    dest='trajectory',
                    help='Trajectory',
                    type=str,
                    required=True)

inputs.add_argument('-sel',
                    action='store',
                    dest='sel',
                    help='Atom selection',
                    type=str,
                    default='protein and resid 262:274')

inputs.add_argument('-dist',
                    action='store',
                    dest='dist',
                    help='Threshold for the cutoff distance',
                    type=str,
                    required=True)

inputs.add_argument('-distNa',
                    action='store',
                    dest='distNa',
                    help='Distance between two Na+ ions',
                    type=str,
                    required=True)

inputs.add_argument('-m',
                    action='store',
                    dest='method',
                    help='bigger, double, inner, outer',
                    type=str,
                    required=True)

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output folder for data',
                    type=str,
                    required=True)

# Parse into useful form
UserInput = parser.parse_args()

u = MDAnalysis.Universe(os.path.join(UserInput.out_dir, UserInput.structure), os.path.join(UserInput.out_dir, UserInput.trajectory))
sodiumLoop=u.select_atoms(UserInput.sel)
sodiumGrp=u.select_atoms('resname SOD')
Inner_site=u.select_atoms('resid 235 and name OD2 OD1 or (resid 236 and name O)')
Outer_site=u.select_atoms('resid 269 272 and name O')
T=len(u.trajectory)
index=[]

# Frames bigger than the input distance
if UserInput.method == 'bigger':
    for i in range(T):
        u.trajectory[i] #Go to current frame
        Na_dist_mean=np.mean(distance.cdist(sodiumLoop.positions,sodiumGrp.positions,'euclidean'),axis=0)
        if np.min(Na_dist_mean)>=float(UserInput.dist):
            index=index + [i+1]            
    np.savetxt(os.path.join(UserInput.out_dir, 'idx' + UserInput.dist + UserInput.method), index, delimiter=' ')

# Frames have two Na+ ions binding
if UserInput.method == 'double':
    for i in range(T):
        u.trajectory[i] #Go to current frame
        Na_dist_mean=np.mean(distance.cdist(sodiumLoop.positions,sodiumGrp.positions,'euclidean'),axis=0)
        Na_min_index=np.argmin(Na_dist_mean)
        Na_min2_index=find_sub_min(Na_dist_mean,2)
        if Na_dist_mean[Na_min2_index]<float(UserInput.dist) and euclidean_distance(sodiumGrp.positions[Na_min_index],sodiumGrp.positions[Na_min2_index]) < float(UserInput.distNa):
            index=index + [i+1]
    np.savetxt(os.path.join(UserInput.out_dir, 'idx' + UserInput.dist + UserInput.method + UserInput.distNa), index, delimiter=' ')

# Frames have inner Na+ ions binding
if UserInput.method == 'inner':
    for i in range(T):
        u.trajectory[i] #Go to current frame
        Na_dist_mean=np.mean(distance.cdist(sodiumLoop.positions,sodiumGrp.positions,'euclidean'),axis=0)
        Inner_dist_mean=np.mean(distance.cdist(Inner_site.positions,sodiumGrp.positions,'euclidean'),axis=0)
        Outer_dist_mean=np.mean(distance.cdist(Outer_site.positions,sodiumGrp.positions,'euclidean'),axis=0)
        Na_min_index=np.argmin(Na_dist_mean)
        Na_min2_index=find_sub_min(Na_dist_mean,2)
        if (np.min(Na_dist_mean)<float(UserInput.dist) and (Na_dist_mean[Na_min2_index]>float(UserInput.dist) or euclidean_distance(sodiumGrp.positions[Na_min_index],sodiumGrp.positions[Na_min2_index]) >= float(UserInput.distNa))) and Inner_dist_mean[Na_min_index] < Outer_dist_mean[Na_min_index]:
            index=index + [i+1]
    np.savetxt(os.path.join(UserInput.out_dir, 'idx' + UserInput.dist + UserInput.method + UserInput.distNa), index, delimiter=' ')

# Frames have one Na+ ions binding
if UserInput.method == 'outer':
    for i in range(T):
        u.trajectory[i] #Go to current frame
        Na_dist_mean=np.mean(distance.cdist(sodiumLoop.positions,sodiumGrp.positions,'euclidean'),axis=0)
        Inner_dist_mean=np.mean(distance.cdist(Inner_site.positions,sodiumGrp.positions,'euclidean'),axis=0)
        Outer_dist_mean=np.mean(distance.cdist(Outer_site.positions,sodiumGrp.positions,'euclidean'),axis=0)
        Na_min_index=np.argmin(Na_dist_mean)
        Na_min2_index=find_sub_min(Na_dist_mean,2)
        if (np.min(Na_dist_mean)<float(UserInput.dist) and (Na_dist_mean[Na_min2_index]>float(UserInput.dist) or euclidean_distance(sodiumGrp.positions[Na_min_index],sodiumGrp.positions[Na_min2_index]) >= float(UserInput.distNa))) and Inner_dist_mean[Na_min_index] > Outer_dist_mean[Na_min_index]:
            index=index + [i+1]
    np.savetxt(os.path.join(UserInput.out_dir, 'idx' + UserInput.dist + UserInput.method + UserInput.distNa), index, delimiter=' ')
