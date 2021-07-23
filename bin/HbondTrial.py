import numpy as np
import argparse
import subprocess
#from shutil import copyfile
import os

#Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description = 'Runs cluster trials over variety of methods and metrics', add_help=False) 

#List all possible user input
inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-s', action='store', dest='structure',help='Structure file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-t', action='store', dest='trajectory',help='Trajectory',type=str,required=True)
inputs.add_argument('-o', action='store', dest='out_dir',help='Output directory',type=str,required=True)

#Parse into useful form
UserInput=parser.parse_args()

# Find Helper scripts
hbond_script = '/home/wud18/python/hbond_analysis_resid.py'

python_cmd = (
    'sbatch --export=t=' + UserInput.trajectory + 
    ',str=' + UserInput.structure +
    ',s=' + hbond_script +
    ',o=' + UserInput.out_dir +
    ' /home/wud18/bash/hbond.slurm'
    )
print(python_cmd)
subprocess.call(python_cmd, shell=True)
