# -*-coding:utf-8-*-
import matplotlib.pyplot as plt
from pylab import text
from PIL import Image
import argparse

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Adding label to images', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-hp', '--help', action='help')
inputs.add_argument('-im',
                    action='store',
                    dest='image',
                    help='Image to add label',
                    type=str,
                    required=True)

inputs.add_argument('-lb',
                    action='store',
                    dest='label',
                    help='label to add',
                    type=str,
                    required=True)

inputs.add_argument('-o',
                    action='store',
                    dest='out_name',
                    help='Output name for the image',
                    type=str,
                    required=True)

# Parse into useful form
UserInput = parser.parse_args()
    
# Loading image
im = Image.open(UserInput.image)
img_size=im.size

# Adding label
fig,ax=plt.subplots(1,figsize=(img_size[0]/100,img_size[1]/100))
ax.imshow(im)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])
ax.text(0.02*img_size[0], 0.06*img_size[1], UserInput.label, size=(img_size[0]+img_size[1])/60, color='black')
fig.savefig(UserInput.out_name)

# Cropping image
im = Image.open(UserInput.out_name)
x=0.11*img_size[0]
y=0.11*img_size[1]
w=0.78*img_size[0]
h=0.78*img_size[1]

# Saving image
region = im.crop((x, y, x+w, y+h))
region.save(UserInput.out_name)
