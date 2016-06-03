#!/usr/bin/env python

import os
import subprocess
import argparse
from PIL import Image
import glob
import signal

# Posix-compliant way to exit shell
signal.signal(signal.SIGTTOU, signal.SIG_IGN)

# Find the helper file generate_shadow.vmd
dir = os.path.dirname(__file__)
cwd = os.getcwd()
shadow_helper = os.path.join(dir, 'generate_shadow.vmd')
middle_helper = os.path.join(dir, 'generate_middle.vmd')
tachyon = os.path.join(dir, 'tachyon')

parser = argparse.ArgumentParser(
        description = (
            'Creates an image with input MIDDLE as solid and SHADOW as shadow.' 
            ), 
        add_help=False
        ) 
#List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-m', '--middle', action='store', dest='middle',help='pdb that should be displayed as solid',type=str,required=True)
inputs.add_argument('-s', '--shadow', action='store', dest='shadow',help='pdb with multiple frames that should be displayed as shadow',type=str,required=True)
inputs.add_argument('-r', '--representation', action='store', dest='rep',help='VMD graphic representation to use',type=str,default='NewCartoon')
inputs.add_argument('-o', '--outdir', action='store', dest='directory',help='output directory',type=str, default=cwd)


#Parse into useful form
UserInput=parser.parse_args()

# Now, let's make some pretty pictures
vmd_render_shadow_cmd = (
        'vmd '
        + UserInput.shadow + ' -dispdev text -e ' 
        + shadow_helper + ' -args '
        + ' -rep ' + UserInput.rep
        )
vmd_render_shadow=subprocess.call([os.getenv('SHELL'), '-i', '-c', vmd_render_shadow_cmd], cwd=UserInput.directory)
os.tcsetpgrp(0,os.getpgrp())

shadow_pattern = UserInput.directory + "/shadow.*.tga"

if os.path.isfile(UserInput.directory + '/shadow.png'):
    os.remove(UserInput.directory + '/shadow.png')

for fname in glob.glob(shadow_pattern):
    shadow_img = Image.open(fname)
    shadow_img = shadow_img.convert("RGBA")
    shadow_data = shadow_img.getdata()

    new_shadow_data = []
    for item in shadow_data:
        if item[0] == 255 and item[1] == 255 and item[2] == 255:
            new_shadow_data.append((255,255,255,0))
        else:
            new_shadow_data.append((item[0], item[1], item[2], 26))

    shadow_img.putdata(new_shadow_data)


    if os.path.isfile(UserInput.directory + '/shadow.png'):
        shadow_old = Image.open(UserInput.directory + '/shadow.png')
        shadow_new = Image.alpha_composite(shadow_img,shadow_old)
        shadow_new.save(UserInput.directory + '/shadow.png')
    else:
        shadow_img.save(UserInput.directory + '/shadow.png',"PNG")

for fname in glob.glob(shadow_pattern):
    os.remove(fname)

vmd_render_middle_cmd = (
        'vmd '
        + UserInput.middle + ' -dispdev text -e ' 
        + middle_helper + ' -args -rep ' + UserInput.rep + ' -outfile '
        + UserInput.directory + '/middle.tga'
        )
vmd_render_middle=subprocess.call([os.getenv('SHELL'), '-i', '-c', vmd_render_middle_cmd], cwd=UserInput.directory)
os.tcsetpgrp(0,os.getpgrp())


#Let's get rid of the white pixels and convert the TGAs to PNGs
middle_img = Image.open(UserInput.directory + '/middle.tga')
middle_img = middle_img.convert("RGBA")
middle_data = middle_img.getdata()

new_middle_data = []
for item in middle_data:
    if item[0] == 89 and item[1] == 89 and item[2] == 89:
        new_middle_data.append((89,89,89,0))
    else:
        new_middle_data.append(item)

middle_img.putdata(new_middle_data)
middle_img.save(UserInput.directory + '/middle.png',"PNG")


#Now, let's layer them together
layered_img = Image.alpha_composite(shadow_new,middle_img)
layered_img.save(UserInput.directory + '/layered.png',"PNG")

blended_img = Image.blend(middle_img,shadow_new,0.5)
blended_img.save(UserInput.directory + '/blended.png',"PNG")



