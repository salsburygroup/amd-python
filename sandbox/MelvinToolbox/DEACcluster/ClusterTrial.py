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
inputs.add_argument('-top', action='store', dest='structure',help='Structure file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-traj', action='store', dest='trajectory',help='Trajectory',type=str,required=True)
inputs.add_argument('-sel', action='store', dest='sel',help='atoms',type=str,required=True)
inputs.add_argument('-o', action='store', dest='out_name',help='Output directory',type=str,required=True)

#Parse into useful form
UserInput=parser.parse_args()

# Find Helper scripts
dir = os.path.dirname(__file__)
cwd = os.getcwd()
cluster_script = '/home/melvrl13/Repos/amd-python/bin/cluster_automatic.py'
python_helpers = dir 


## AP
if not os.path.exists(os.path.join(UserInput.out_name, 'HDBSCAN')):
    os.makedirs(os.path.join(UserInput.out_name, 'HDBSCAN'))
python_HDBSCAN_cmd = (
    'sbatch --export=t=' + UserInput.trajectory + ',str=' + UserInput.structure +
    ',s=' + cluster_script +
    ',m=HDBSCAN' + 
    ',a=' + UserInput.sel + 
    ',o=' + os.path.join(UserInput.out_name, 'HDBSCAN') + 
    ' ' +  os.path.join(python_helpers, 'PythonSubmit.slurm')
    )
print(python_HDBSCAN_cmd)
subprocess.call(python_HDBSCAN_cmd, shell=True)


## AP
if not os.path.exists(os.path.join(UserInput.out_name, 'IMWKRescaled')):
    os.makedirs(os.path.join(UserInput.out_name, 'IMWKRescaled'))
python_IMWKRescaled_cmd = (
    'sbatch --export=t=' + UserInput.trajectory + ',str=' + UserInput.structure +
    ',s=' + cluster_script +
    ',m=IMWKRescaled' + 
    ',a=' + UserInput.sel + 
    ',o=' + os.path.join(UserInput.out_name, 'IMWKRescaled') + 
    ' ' +  os.path.join(python_helpers, 'PythonSubmit.slurm')
    )
print(python_IMWKRescaled_cmd)
subprocess.call(python_IMWKRescaled_cmd, shell=True)

## AP
#if not os.path.exists(os.path.join(UserInput.out_name, 'AffinityPropagation')):
#    os.makedirs(os.path.join(UserInput.out_name, 'AffinityPropagation'))
#python_AffinityPropagation_cmd = (
#    'sbatch --export=t=' + UserInput.trajectory + ',str=' + UserInput.structure +
#    ',s=' + cluster_script +
#    ',m=AffinityPropagation' + 
#    ',a=' + UserInput.sel + 
#    ',o=' + os.path.join(UserInput.out_name, 'AffinityPropagation') + 
#    ' ' +  os.path.join(python_helpers, 'PythonSubmit.slurm')
#    )
#print(python_AffinityPropagation_cmd)
#subprocess.call(python_AffinityPropagation_cmd, shell=True)
#
### GMM
#if not os.path.exists(os.path.join(UserInput.out_name, 'GMM')):
#    os.makedirs(os.path.join(UserInput.out_name, 'GMM'))
#python_GMM_cmd = (
#    'sbatch --export=t=' + UserInput.trajectory + ',str=' + UserInput.structure +
#    ',s=' + cluster_script +
#    ',m=GMM'
#    ',a=' + UserInput.sel + 
#    ',o=' + os.path.join(UserInput.out_name, 'GMM') + 
#    ' ' +  os.path.join(python_helpers, 'PythonSubmit.slurm')
#    )
#print(python_GMM_cmd)
#subprocess.call(python_GMM_cmd, shell=True)
#
### KMeans
#if not os.path.exists(os.path.join(UserInput.out_name, 'KMeans')):
#    os.makedirs(os.path.join(UserInput.out_name, 'KMeans'))
#python_KMeans_cmd = (
#    'sbatch --export=t=' + UserInput.trajectory + ',str=' + UserInput.structure +
#    ',s=' + cluster_script +
#    ',m=KMeans'
#    ',a=' + UserInput.sel + 
#    ',o=' + os.path.join(UserInput.out_name, 'KMeans') + 
#    ' ' +  os.path.join(python_helpers, 'PythonSubmit.slurm')
#    )
#print(python_KMeans_cmd)
#subprocess.call(python_KMeans_cmd, shell=True)
#
### MS
#if not os.path.exists(os.path.join(UserInput.out_name, 'MeanShift')):
#    os.makedirs(os.path.join(UserInput.out_name, 'MeanShift'))
#python_MeanShift_cmd = (
#    'sbatch --export=t=' + UserInput.trajectory + ',str=' + UserInput.structure +
#    ',s=' + cluster_script +
#    ',m=MeanShift'
#    ',a=' + UserInput.sel + 
#    ',o=' + os.path.join(UserInput.out_name, 'MeanShift') + 
#    ' ' +  os.path.join(python_helpers, 'PythonSubmit.slurm')
#    )
#print(python_MeanShift_cmd)
#subprocess.call(python_MeanShift_cmd, shell=True)
#
### VBGMM
#if not os.path.exists(os.path.join(UserInput.out_name, 'VBGMM')):
#    os.makedirs(os.path.join(UserInput.out_name, 'VBGMM'))
#python_VBGMM_cmd = (
#    'sbatch --export=t=' + UserInput.trajectory + ',str=' + UserInput.structure +
#    ',s=' + cluster_script +
#    ',m=VBGMM'
#    ',a=' + UserInput.sel + 
#    ',o=' + os.path.join(UserInput.out_name, 'VBGMM') + 
#    ' ' +  os.path.join(python_helpers, 'PythonSubmit.slurm')
#    )
#print(python_VBGMM_cmd)
#subprocess.call(python_VBGMM_cmd, shell=True)
