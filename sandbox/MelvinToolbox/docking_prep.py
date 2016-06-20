def replace_atoms(input,template,output):

    import re

    # Open file with only #ATOM lines
    initial = open(input,'r+').read()
    initial = initial.rstrip()

    # Create dictionary 
    key=re.findall('(ATOM\s*[0-9]+).*\n',initial)
    value=re.findall('(ATOM\s*[0-9]+.*)\n',initial)
    replacements = dict(zip(key,value))

    with open(template,'rb') as file, open(output,'w') as out:
        for line in file:
            test = re.findall('(ATOM\s*[0-9]+).*\n',line)
            if test:
                text = str(replacements[test[0]]) +'\n'
            else: 
                text = line
            out.write(text)

