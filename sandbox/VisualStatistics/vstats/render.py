def shadow(rep, number, directory):

	#Get path to vmd helper
	dir = os.path.dirname(__file__)
	shadow_helper = os.path.join(dir, 'vmd/generate_shadow.vmd')

	#Render the shadow
	vmd_render_shadow_cmd = (os.environ["VMD_HOME"] + '/vmd_MACOSXX86 ' 
		+ directory +'/sigma.pdb -dispdev text -e ' 
		+  + ' -args -first 1 -last ' + str(len(number)) + ' -rep ' + rep + ' -outfile ' 
		+ directory + '/shadow.dat')
	vmd_render_shadow=subprocess.call(vmd_render_shadow_cmd,shell=True)

	tachyon_render_shadow_cmd = (os.environ["VMD_HOME"] + '/tachyon_MACOSXX86 ' +
		' -trans_vmd ' + directory + '/shadow.dat -o ' + directory + '/shadow.tga')
	tachyon_render_shadow=subprocess.call(tachyon_render_shadow_cmd,shell=True)

	#Now, let's work on the alpha channel
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

def solid(rep, directory):
	
	#Get path to vmd helper
	dir = os.path.dirname(__file__)
	solid_helper = os.path.join(dir, 'vmd/generate_solid.vmd')

	#Render the center/ solid image
	vmd_render_solid_cmd = (os.environ["VMD_HOME"] + '/vmd_MACOSXX86 ' 
		+ directory +'/solid.pdb -dispdev text -e ' 
		+ solid_helper + ' -args -outfile '
		+directory + '/solid.dat' + ' -rep ' + rep)
	vmd_render_solid=subprocess.call(vmd_render_solid_cmd,shell=True)

	tachyon_render_solid_cmd = (os.environ["VMD_HOME"] + '/tachyon_MACOSXX86'+
                ' -trans_vmd ' + directory + '/solid.dat -o ' + directory + '/solid.tga')
	tachyon_render_solid=subprocess.call(tachyon_render_solid_cmd,shell=True)

	#Now, let's work on the alpha channel
	median_img.putdata(new_median_data)
        median_img.save(directory + '/median.png',"PNG")
	new_median_data = []
        for item in median_data:
		if item[0] == 89 and item[1] == 89 and item[2] == 89:
			new_median_data.append((89,89,89,0))
		else:
			new_median_data.append(item)

	median_img.putdata(new_median_data)
	median_img.save(directory + '/median.png',"PNG")

def layer(img1, img2, directory):

	layered_img = Image.alpha_composite(img2,img1)
	layered_img.save(directory + '/layered.png',"PNG")

def blend(img1, img2, directory):

	blended_img = Image.blend(img1, img2, 0.5)
	blended_img.save(directory + '/blended.png',"PNG")
