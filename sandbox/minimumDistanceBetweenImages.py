# coding: utf-8
#get_ipython().magic(u'matplotlib inline')
#!/usr/bin/env python
import MDAnalysis
import numpy as np
import MDAnalysis.analysis.distances
import argparse
import matplotlib.pyplot as plt

# Jiajie Xiao
# Oct 2, 2016

parser = argparse.ArgumentParser(
    description='Quantification of the geometric relationship between HIS72 and THY3',
    add_help=False
)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')

inputs.add_argument(
    '-structure',
    action='store',
    dest='structure',
    help='Structure file corresponding to trajectory',
    type=str,
    required=True
)

inputs.add_argument(
    '-traj',
    action='store',
    dest='trajectory',
    help='Trajectory',
    type=str,
    required=True
)

inputs.add_argument(
    '-stride',
    action='store',
    dest='stride',
    help='stride',
    type=int,
    default=1
)

#inputs.add_argument(
#    '-sel',
#    action='store',
#    dest='sel',
#    help='Atom selection',
#    type=str,
#    default='all'
#)

#inputs.add_argument(
#    '-o',
#    action='store',
#    dest='out_name',
#    help='Output prefix',
#    type=str,
#    required=True
#)

# Parse into useful form
UserInput = parser.parse_args()

u = MDAnalysis.Universe(UserInput.structure, UserInput.trajectory)


# In[67]:

#u.dimensions


# In[73]:

protein=u.select_atoms('protein')


# In[74]:

nucleic=u.select_atoms('nucleic')


# In[77]:

TBAcomplex=u.select_atoms('nucleic or protein')


# In[75]:

#protein.center_of_mass()


# In[76]:

#nucleic.center_of_mass()


# In[78]:

#u.atoms.center_of_geometry()

T=len(u.trajectory)
#T=5000
data = np.zeros((T/UserInput.stride,4))
j=0

for i in range(0,T,UserInput.stride):
    u.trajectory[i]
    d0=np.linalg.norm(u.atoms.center_of_geometry()-TBAcomplex.center_of_mass()) # distance between the center of system and COM of the TBAcomplex
    d1=np.min(MDAnalysis.analysis.distances.distance_array((TBAcomplex.positions+[u.dimensions[0],0,0]).astype('float32'),TBAcomplex.positions)) 
    d2=np.min(MDAnalysis.analysis.distances.distance_array((TBAcomplex.positions+[0,u.dimensions[1],0]).astype('float32'),TBAcomplex.positions))
    d3=np.min(MDAnalysis.analysis.distances.distance_array((TBAcomplex.positions+[0,0,u.dimensions[2]]).astype('float32'),TBAcomplex.positions))
    data[j] = [d0,d1,d2,d3]
    j=j+1

plt.figure()
plt.plot(data[:,0]) #plot only the second column
plt.xlabel('Frame')
plt.ylabel(r'Distance($\AA$)')
plt.title('Distance between the COM of TBAcomplex and geometry center of simulation box')
#plt.style.use('bmh') #Uncomment to use matplotlib style sheets. Also, remove "import seaborn."
plt.savefig('minimumDistanceBetweenImages_1.png',transparent=True)
plt.savefig('minimumDistanceBetweenImages_1.pdf',transparent=True)

plt.figure()
plt.plot(data[:,1:4]) #plot only the second column
plt.xlabel('Frame')
plt.ylabel(r'Distance($\AA$)')
plt.title('Shortest distance between TBAcomplex and its image')
#plt.style.use('bmh') #Uncomment to use matplotlib style sheets. Also, remove "import seaborn."
plt.savefig('minimumDistanceBetweenImages_2.png',transparent=True)
plt.savefig('minimumDistanceBetweenImages_2.pdf',transparent=True)

np.savetxt('minimumDistanceBetweenImages.dat',data)
# In[46]:

#protein.bbox()


# In[79]:

#nucleic.bbox()


# In[80]:



# In[81]:



# In[82]:



# In[83]:



# In[84]:

#(TBAcomplex.positions+[u.dimensions[0],0,0]).dtype


# In[85]:

# u.atoms.pack_into_box 


# In[86]:

#u.atoms.write("TBAcomplex_K_565frames_stride100_packintobox.pdb")


# In[ ]:



