import numpy as np
import argparse
import subprocess
from shutil import copyfile
import os

#Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description = 'Runs cluster trials over variety of methods and metrics', add_help=False) 

#List all possible user input
inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-p', action='store', dest='path',help='Path of the directory',type=str,required=True)
inputs.add_argument('-s', action='store', dest='structure',help='Structure file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-t', action='store', dest='trajectory',help='Trajectory',type=str,required=True)
inputs.add_argument('-title', action='store', dest='title',help='Title',type=str,required=True)
inputs.add_argument('-o', action='store', dest='out_name',help='Output directory',type=str,required=True)

# Parse into useful form
UserInput=parser.parse_args()

## SASA stride 1
if not os.path.exists(os.path.join(UserInput.out_name, 'whole_trajectory')):
    os.makedirs(os.path.join(UserInput.out_name, 'whole_trajectory'))
python_SASA_cmd = (
    'sbatch --export=t=' + UserInput.trajectory + 
    ',str=' + UserInput.structure +
    ',b=' + UserInput.title + 
    ',o=' + os.path.join(UserInput.out_name, 'whole_trajectory') + 
    ' /home/wud18/bash/SASA.slurm'
    )
print(python_SASA_cmd)
subprocess.call(python_SASA_cmd, shell=True)

## SASA stride 10
if not os.path.exists(os.path.join(UserInput.out_name, 'trajectory_stride10')):
    os.makedirs(os.path.join(UserInput.out_name, 'trajectory_stride10'))
python_SASA_stride10_cmd = (
    'sbatch --export=t=' + UserInput.path + '/protein_stride10_aligned.dcd' +
    ',str=' + UserInput.structure +
    ',b=' + UserInput.title + 
    ',o=' + os.path.join(UserInput.out_name, 'trajectory_stride10') + 
    ' /home/wud18/bash/SASA_stride10.slurm'
    )
print(python_SASA_stride10_cmd)
subprocess.call(python_SASA_stride10_cmd, shell=True)
