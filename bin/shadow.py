import time
start = time.clock()
import os
import PIL.Image
import glob
import argparse

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Make shadow.png', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output directory',
                    type=str,
                    required=True)

# Parse into useful form
UserInput = parser.parse_args()

shadow_pattern = UserInput.out_dir + "/shadow.*.tga"
if os.path.isfile(UserInput.out_dir + '/shadow.png'):
    os.remove(UserInput.out_dir + '/shadow.png')
    
for file in glob.glob(shadow_pattern):
    shadow_img = PIL.Image.open(file)
    shadow_img = shadow_img.convert("RGBA")
    shadow_data = shadow_img.getdata()
        
    new_shadow_data = []
    for item in shadow_data:
        if item[0] == 255 and item[1] == 255 and item[2] == 255:
            new_shadow_data.append((255, 255, 255, 0))
        else:
            new_shadow_data.append((item[0], item[1], item[2], 26))
                
    shadow_img.putdata(new_shadow_data)
            
    if os.path.isfile(UserInput.out_dir + '/shadow.png'):
        shadow_old = PIL.Image.open(UserInput.out_dir + '/shadow.png')
        shadow_new = PIL.Image.alpha_composite(shadow_img, shadow_old)
        shadow_new.save(UserInput.out_dir + '/shadow.png')
    else:
        shadow_img.save(UserInput.out_dir + '/shadow.png', "PNG")
            
#for file in glob.glob(shadow_pattern):
#    os.remove(file)
end = time.clock()
print(end - start)
