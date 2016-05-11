#!/usr/env Python

#Ryan Melvin
#Tue Jun 30 10:56:11 EDT 2015

#Outputs plot with same base name as RMSD file in same path as RMSD dat file

import numpy as np
import matplotlib.pyplot as plt
import argparse
import os.path as path
import seaborn

parser = argparse.ArgumentParser(description = 'Plot column 2 of VMD-generated trajrmsd.dat file', add_help=False) 

inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-trajrmsd', action='store', dest='rmsd_file',help='VMD-generated RMSD text file path as string',type=str,required=True)
inputs.add_argument('-title', action='store', dest='plt_title', help='Title for plot', type=str, default='RMSD time series')

UserInput=parser.parse_args()

# Get working directory
wd = path.dirname(UserInput.rmsd_file)

# Get file name 
base = path.basename(UserInput.rmsd_file)
outname = path.splitext(base)[0]

# Load text file. Skip headers and reference frame 1
rmsd = np.loadtxt(UserInput.rmsd_file,skiprows=2)

#plot
plt.figure()
plt.plot(rmsd[:,1]) #plot only the second column
plt.xlabel('Frame')
plt.ylabel('RMSD')
plt.title(UserInput.plt_title)
#plt.style.use('bmh') #Uncomment to use matplotlib style sheets. Also, remove "import seaborn."

plt.savefig(wd + '/' + outname + '.png')
