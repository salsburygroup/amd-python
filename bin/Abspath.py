import os

name=os.popen("pwd")
print(name.readlines()[-1].split('\n')[0])
