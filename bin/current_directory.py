import os

#name=os.popen("pwd")
#print(name.readlines()[0].strip().split('/')[-1])
print(os.getcwd())
print(os.getcwd().split('/')[-1])
