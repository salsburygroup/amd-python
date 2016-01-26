#!/usr/bin/env python 

# Titrates a molecule with the Zinc ion

# Import dependencies
import subprocess
import os
import shutil
import glob

# Parameters
Device = 2
Max = 5

# Paths on GPU
analysis_scripts = '/opt/AMD/amd-vmd/VMDScripts/Analysis/'
setup_scripts = '/opt/AMD/amd-vmd/VMDScripts/Setup/'
vmd_exec = '/usr/local/bin/vmd'

# Locate helper files
dir = os.path.dirname(__file__)
last_frame_helper = os.path.join(dir,'last_frame.vmd')
ionize_helper = os.path.join(dir,'ionize.vmd')

# Locate acemd parameters
param = glob.glob('*.prm')[0]
    
# Get working directory 
cwd = os.getcwd()

# First run

for i in range(0,Max+1):
    dir_to_make = os.path.join(cwd,str(i))
    if not os.path.exists(dir_to_make):
        os.makedirs(dir_to_make)
    # Move files
    shutil.copy(os.path.join(cwd,'structure.psf'),os.path.join(dir_to_make,'structure.psf'))
    shutil.copy(os.path.join(cwd,'structure.pdb'),os.path.join(dir_to_make,'structure.pdb'))
    shutil.copy(os.path.join(cwd,param),os.path.join(dir_to_make,param))
    shutil.copy(os.path.join(cwd,'Titrate_Relax.conf'),os.path.join(dir_to_make,'Titrate_Relax.conf'))

    #Run simulation in new folder
    os.chdir(dir_to_make)
    acemd_cmd = 'acemd --device ' + str(Device) +' Titrate_Relax.conf >&relax.log' 
    process = subprocess.Popen(acemd_cmd, shell=True)
    process.wait()
    
    # Wrap a trajectory
    vmd_wrap_cmd = (
    vmd_exec + ' structure.psf out.dcd -dispdev text -e ' +  
        analysis_scripts + 'Wrap.tcl -args -atomsel ' 
        + 'segid_F10_and_noh' + ' -outfile out.dcd'
            ) 
    vmd_wrap = subprocess.call(vmd_wrap_cmd, shell=True)
    
    # Get last frame of trajectory
    vmd_last_frame_cmd = (
        vmd_exec + 
        ' structure.psf out.dcd -dispdev text -e ' +  last_frame_helper
            ) 
    vmd_last_frame = subprocess.call(vmd_last_frame_cmd, shell=True)

    # Add next ion in titration
    vmd_ionize_cmd = (vmd_exec + ' -dispdev text -e ' + ionize_helper)
    vmd_ionize = subprocess.call(vmd_ionize_cmd,shell=True)

    shutil.copy('next.psf',os.path.join(cwd,'structure.psf'))
    shutil.copy('next.pdb',os.path.join(cwd,'structure.pdb'))

