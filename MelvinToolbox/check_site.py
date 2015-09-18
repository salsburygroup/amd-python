#!/usr/bin/env python

#See if F10+acid is in expected acid binding site on Albumin

import MDAnalysis as md
import argparse
import csv

# Initialize argument parser 
parser = argparse.ArgumentParser()

# Get PDB to check from user input
parser.add_argument('-pdb', action='store', dest='p',type=str,required=True)

# Get path to output file
parser.add_argument('-o', action='store', dest='o',type=str,required=True)

inargs = parser.parse_args()
pdb_str = inargs.p
outfile = inargs.o

# Load F10-acid for comparison
F10_acid = md.Universe(pdb_str)

#site1 = albumin.selectAtoms("resnum 1001") 
site1 = '34.13238144 14.35730648 36.12953949'
site2 = '49.07566833 12.52499866 18.63955498'
site3 = '11.25188351 11.23558807 15.06564522'
site4 = '10.6989994 2.22866702 20.79711151'
site5 = '-0.7801764 6.19541216 39.39217758' 
site6 = '23.67372322 3.62183332 -1.2207222'
site7 = '33.92470169 15.60519981 9.01439857'


models_in_site = []

for frame in range(len(F10_acid.trajectory)):
    F10_acid.trajectory[frame]
    in_site = F10_acid.selectAtoms('(point %s 5) or (point %s 5) or (point %s 5) or (point %s 5) or (point %s 5) or (point %s 5) or (point %s 5)' % (site1, site2, site3, site4, site5, site6, site7))
    if len(in_site)>0:
        model_number = frame + 1
        models_in_site.append(model_number)

# Write list of models to output file
with open(outfile, 'wb') as file:
    wr = csv.writer(file)
    wr.writerow(models_in_site)
