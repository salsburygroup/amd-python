#/usr/bin/env python
from PIL import Image
#Let's get rid of the white pixels and convert the TGAs to PNGs
middle_img = Image.open('middle.tga')
middle_img = middle_img.convert("RGBA")
middle_data = middle_img.getdata()

new_middle_data = []
for item in middle_data:
    if item[0] == 255 and item[1] == 255 and item[2] == 255:
        new_middle_data.append((255,255,255,0))
    else:
        new_middle_data.append(item)

middle_img.putdata(new_middle_data)
middle_img.save('middle.png',"PNG")

shadow_img = Image.open('shadow.tga')
shadow_img = shadow_img.convert("RGBA")
shadow_data = shadow_img.getdata()

new_shadow_data = []
for item in shadow_data:
    if item[0] == 255 and item[1] == 255 and item[2] == 255:
        new_shadow_data.append((255,255,255,0))
    else:
        new_shadow_data.append(item)

shadow_img.putdata(new_shadow_data)
shadow_img.save('shadow.png',"PNG")

#Now, let's layer them together
layered_img = Image.alpha_composite(shadow_img,middle_img)
layered_img.save('layered.png',"PNG")

blended_img = Image.blend(middle_img,shadow_img,0.7)
blended_img.save('blended.png',"PNG")

