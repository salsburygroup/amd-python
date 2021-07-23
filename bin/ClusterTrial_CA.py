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
inputs.add_argument('-s', action='store', dest='structure',help='Structure file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-t', action='store', dest='trajectory',help='Trajectory',type=str,required=True)
inputs.add_argument('-title', action='store', dest='title',help='Trajectory',type=str,required=True)
inputs.add_argument('-sel', action='store', dest='sel',help='atoms',type=str,required=True)
inputs.add_argument('-o', action='store', dest='out_name',help='Output directory',type=str,required=True)
inputs.add_argument('-tm', action='store', dest='timestep', help='Timestep between two frames', type=str, required=True)

#Parse into useful form
UserInput=parser.parse_args()

# Find Helper scripts
cwd = os.getcwd()
cluster_script = '/home/wud18/python/cluster_automatic.py'

## HDBSCAN
if not os.path.exists(os.path.join(UserInput.out_name, 'HDBSCAN')):
    os.makedirs(os.path.join(UserInput.out_name, 'HDBSCAN'))
python_HDBSCAN_cmd = (
    'sbatch --export=t=' + UserInput.trajectory + ',str=' + UserInput.structure +
    ',s=' + cluster_script +
    ',m=HDBSCAN' +
    ',a=' + UserInput.sel + 
    ',b=' + UserInput.title + 
    ',tm=' + UserInput.timestep + 
    ',o=' + os.path.join(UserInput.out_name, 'HDBSCAN') + 
    ' /home/wud18/bash/HDsubmit_CA.slurm'
    )
print(python_HDBSCAN_cmd)
subprocess.call(python_HDBSCAN_cmd, shell=True)

## IMWKRescaled
if not os.path.exists(os.path.join(UserInput.out_name, 'IMWKRescaled')):
    os.makedirs(os.path.join(UserInput.out_name, 'IMWKRescaled'))
python_IMWKRescaled_cmd = (
    'sbatch --export=t=' + UserInput.trajectory + ',str=' + UserInput.structure +
    ',s=' + cluster_script +
    ',m=IMWKRescaled' + 
    ',a=' + UserInput.sel + 
    ',b=' + UserInput.title + 
    ',tm=' + UserInput.timestep + 
    ',o=' + os.path.join(UserInput.out_name, 'IMWKRescaled') + 
    ' /home/wud18/bash/AHsubmit_CA.slurm'
    )
print(python_IMWKRescaled_cmd)
subprocess.call(python_IMWKRescaled_cmd, shell=True)
