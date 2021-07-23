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

# Parse into useful form
UserInput = parser.parse_args()
    
# Loading image
im = Image.open(UserInput.image)

# Image size
img_size = im.size
print("Image size:{}".format(img_size))

