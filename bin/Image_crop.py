# -*-coding:utf-8-*-
from PIL import Image
import argparse

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Cropping images', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-hp', '--help', action='help')
inputs.add_argument('-im',
                    action='store',
                    dest='image',
                    help='Image to crop',
                    type=str,
                    required=True)

inputs.add_argument('-x',
                    action='store',
                    dest='x',
                    help='Pixels in x-axis',
                    type=str,
                    required=True)

inputs.add_argument('-y',
                    action='store',
                    dest='y',
                    help='Pixels in y-axis',
                    type=str,
                    required=True)

inputs.add_argument('-w',
                    action='store',
                    dest='width',
                    help='Image width in pixels',
                    type=str,
                    required=True)

inputs.add_argument('-h',
                    action='store',
                    dest='height',
                    help='Image height in pixels',
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

# Cropping image
x = int(UserInput.x)
y = int(UserInput.y)
w = int(UserInput.width)
h = int(UserInput.height)

# Saving image
region = im.crop((x, y, x+w, y+h))
region.save(UserInput.out_name)

