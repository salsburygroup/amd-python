#! /usr/bin/env python
# Visual distribution generator

def vdistribution(middle, shadow, number, rep, directory):
    import os
    import subprocess
    from PIL import Image
    from sys import exit

    # Find the helper file generate_shadow.vmd
    dir = os.path.dirname(__file__)
    shadow_helper = os.path.join(dir, 'generate_shadow.vmd')
    middle_helper = os.path.join(dir, 'generate_middle.vmd')

    #Get path to VMD
    if "VMD_HOME" in os.environ:
        vmd_path = os.environ["VMD_HOME"]
    else:
        exit("VMD_HOME not found in environmental variables.")

    # Now, let's make some pretty pictures
    vmd_render_shadow_cmd = (
            vmd_path + '/vmd_MACOSXX86 ' 
            + shadow + ' -dispdev text -e ' 
            + shadow_helper + ' -args -first 1 -last ' + str(number -1 ) 
            + ' -rep ' + rep + ' -outfile ' + directory + '/shadow.dat'
            )
    vmd_render_shadow=subprocess.call(vmd_render_shadow_cmd,shell=True)


    vmd_render_middle_cmd = (
            vmd_path + '/vmd_MACOSXX86 ' 
            + middle + ' -dispdev text -e ' 
            + middle_helper + ' -args -rep ' + rep + ' -outfile '
            + directory + '/middle.dat'
            )
    vmd_render_middle=subprocess.call(vmd_render_middle_cmd,shell=True)

    tachyon_render_shadow_cmd = (vmd_path + '/tachyon_MACOSXX86'+
            ' -trans_vmd ' + directory + '/shadow.dat -o ' + directory + '/shadow.tga')
    tachyon_render_shadow=subprocess.call(tachyon_render_shadow_cmd,shell=True)

    tachyon_render_middle_cmd = (vmd_path + '/tachyon_MACOSXX86'+
            ' -trans_vmd ' + directory + '/middle.dat -o ' + directory + '/middle.tga')
    tachyon_render_middle=subprocess.call(tachyon_render_middle_cmd,shell=True)

    #Let's get rid of the white pixels and convert the TGAs to PNGs
    middle_img = Image.open(directory + '/middle.tga')
    middle_img = middle_img.convert("RGBA")
    middle_data = middle_img.getdata()

    new_middle_data = []
    for item in middle_data:
        if item[0] == 89 and item[1] == 89 and item[2] == 89:
            new_middle_data.append((89,89,89,0))
        else:
            new_middle_data.append(item)

    middle_img.putdata(new_middle_data)
    middle_img.save(directory + '/middle.png',"PNG")

    shadow_img = Image.open(directory + '/shadow.tga')
    shadow_img = shadow_img.convert("RGBA")
    shadow_data = shadow_img.getdata()

    new_shadow_data = []
    for item in shadow_data:
        if item[0] == 255 and item[1] == 255 and item[2] == 255:
            new_shadow_data.append((255,255,255,0))
        else:
            new_shadow_data.append(item)

    shadow_img.putdata(new_shadow_data)
    shadow_img.save(directory + '/shadow.png',"PNG")

    #Now, let's layer them together
    layered_img = Image.alpha_composite(shadow_img,middle_img)
    layered_img.save(directory + '/layered.png',"PNG")

    blended_img = Image.blend(middle_img,shadow_img,0.5)
    blended_img.save(directory + '/blended.png',"PNG")

if __name__ == '__main__':
    import os
    import argparse
    cwd = os.getcwd()

    parser = argparse.ArgumentParser(
            description = (
                'Outputs layered and blended image of the specified conformations' 
                ), 
            add_help=False
            ) 
    #List all possible user input
    inputs = parser.add_argument_group('Input arguments')
    inputs.add_argument('-h', '--help', action='help')
    inputs.add_argument('-m', '--middle', action='store', dest='middle',help='pdb that should be displayed as solid',type=str,required=True)
    inputs.add_argument('-s', '--shadow', action='store', dest='shadow',help='pdb with multiple frames that should be displayed as shadow',type=str,required=True)
    inputs.add_argument('-n', '--number', action='store', dest='number',help='number of frames in shadow',type=int,required=True)
    inputs.add_argument('-r', '--representation', action='store', dest='rep',help='VMD graphic representation to use',type=str,default='NewRibbons')
    inputs.add_argument('-o', '--outdir', action='store', dest='directory',help='output directory',type=str, default=cwd)


    #Parse into useful form
    UserInput=parser.parse_args()

    vdistribution(UserInput.middle, UserInput.shadow, UserInput.number, UserInput.rep, UserInput.directory)



