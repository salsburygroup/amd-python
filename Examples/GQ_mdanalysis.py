# -*- coding: utf-8 -*-
#!/usr/bin/env python
import numpy as np
import MDAnalysis
import argparse
import matplotlib.pyplot as plt

# Jiajie Xiao, Aug 8th 2016
# CITATION DOI:10.1080/07391102.2015.1055303

parser = argparse.ArgumentParser(
    description='Return quantities systematically describing conformational rearragements' +
                '  in G-quadruplexes.',
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

#Quartet numeration begins at the 5′-terminus. The guani- nes in the 1st quartet are numbered in the direction of the hydrogen bonding polarity, from the donor to the accep- tor, and from the 5′-terminus. From http://www.bmsc.washington.edu/CrystaLinks/man/pdb/part_62.html Nucleic acid residues are listed from the 5' to the 3' terminus.



# Upper G-quartet
GQ_1_1 = u.select_atoms("resid 1 and nucleicbase")
GQ_1_2 = u.select_atoms("resid 6 and nucleicbase")
GQ_1_3 = u.select_atoms("resid 10 and nucleicbase")
GQ_1_4 = u.select_atoms("resid 15 and nucleicbase")

# Lower  G-quartet
GQ_2_1 = u.select_atoms("resid 2 and nucleicbase")
GQ_2_2 = u.select_atoms("resid 5 and nucleicbase")
GQ_2_3 = u.select_atoms("resid 11 and nucleicbase")
GQ_2_4 = u.select_atoms("resid 14 and nucleicbase")



T=len(u.trajectory)
#T=10
dist_GQ1_GQ2 = np.zeros((T,1))
dist_GQ_Gua = np.zeros((T, 8))
twistAngle = np.zeros((T,1))
dihedral = np.zeros((T, 8))
phi_quartet_axis_norm_Gua = np.zeros((T, 8))
omega = np.zeros((T, 12))

for i in range(T):
    u.trajectory[i]

    # Distance between upper and lower quartets
    dist_GQ1_GQ2[i] = np.linalg.norm((GQ_1_1+GQ_1_2+GQ_1_3+GQ_1_4).center_of_mass()
    -(GQ_2_1+GQ_2_2+GQ_2_3+GQ_2_4).center_of_mass())
               
    # The distances between the quartet COM and COMs of particular guanines
    dist_GQ1_GQ_1_1 = np.linalg.norm((GQ_1_1+GQ_1_2+GQ_1_3+GQ_1_4).center_of_mass()
    -GQ_1_1.center_of_mass())
    dist_GQ1_GQ_1_2 = np.linalg.norm((GQ_1_1+GQ_1_2+GQ_1_3+GQ_1_4).center_of_mass()
    -GQ_1_2.center_of_mass())
    dist_GQ1_GQ_1_3 = np.linalg.norm((GQ_1_1+GQ_1_2+GQ_1_3+GQ_1_4).center_of_mass()
    -GQ_1_3.center_of_mass())
    dist_GQ1_GQ_1_4 = np.linalg.norm((GQ_1_1+GQ_1_2+GQ_1_3+GQ_1_4).center_of_mass()
    -GQ_1_4.center_of_mass())

    dist_GQ2_GQ_2_1 = np.linalg.norm((GQ_2_1+GQ_2_2+GQ_2_3+GQ_2_4).center_of_mass()
    -GQ_2_1.center_of_mass())
    dist_GQ2_GQ_2_2 = np.linalg.norm((GQ_2_1+GQ_2_2+GQ_2_3+GQ_2_4).center_of_mass()
    -GQ_2_2.center_of_mass())
    dist_GQ2_GQ_2_3 = np.linalg.norm((GQ_2_1+GQ_2_2+GQ_2_3+GQ_2_4).center_of_mass()
    -GQ_2_3.center_of_mass())
    dist_GQ2_GQ_2_4 = np.linalg.norm((GQ_2_1+GQ_2_2+GQ_2_3+GQ_2_4).center_of_mass()
    -GQ_2_4.center_of_mass())

    dist_GQ_Gua [i] = [dist_GQ1_GQ_1_1, dist_GQ1_GQ_1_2, dist_GQ1_GQ_1_3, dist_GQ1_GQ_1_4,
    dist_GQ2_GQ_2_1, dist_GQ2_GQ_2_2, dist_GQ2_GQ_2_3, dist_GQ2_GQ_2_4]
               
    # Twist angle between quatet (polarity) 
    v1 = (GQ_1_1.select_atoms("name N9").coordinates()
    -GQ_1_2.select_atoms("name N9").coordinates())[0]
    v2 = (GQ_2_1.select_atoms("name N9").coordinates()
    -GQ_2_2.select_atoms("name N9").coordinates())[0]
    twistAngle_1=180/np.pi*np.arccos(np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))

    v1 = (GQ_1_2.select_atoms("name N9").coordinates()
    -GQ_1_3.select_atoms("name N9").coordinates())[0]
    v2 = (GQ_2_2.select_atoms("name N9").coordinates()
    -GQ_2_3.select_atoms("name N9").coordinates())[0]
    twistAngle_2=180/np.pi*np.arccos(np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))

    v1 = (GQ_1_3.select_atoms("name N9").coordinates()
    -GQ_1_4.select_atoms("name N9").coordinates())[0]
    v2 = (GQ_2_3.select_atoms("name N9").coordinates()
    -GQ_2_4.select_atoms("name N9").coordinates())[0]
    twistAngle_3=180/np.pi*np.arccos(np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))

    v1 = (GQ_1_4.select_atoms("name N9").coordinates()
    -GQ_1_1.select_atoms("name N9").coordinates())[0]
    v2 = (GQ_2_4.select_atoms("name N9").coordinates()
    -GQ_2_1.select_atoms("name N9").coordinates())[0]
    twistAngle_4=180/np.pi*np.arccos(np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))

    twistAngle[i] = (twistAngle_1+twistAngle_2+twistAngle_3+twistAngle_4)/4

    dihedral_1_1234 = MDAnalysis.core.topologyobjects.Dihedral(
    GQ_1_1.select_atoms("name N1")+GQ_1_2.select_atoms("name N1")
    +GQ_1_3.select_atoms("name N1")+GQ_1_4.select_atoms("name N1")).value()

    dihedral_1_2341 = MDAnalysis.core.topologyobjects.Dihedral(
    GQ_1_2.select_atoms("name N1")+GQ_1_3.select_atoms("name N1")
    +GQ_1_4.select_atoms("name N1")+GQ_1_1.select_atoms("name N1")).value()

    dihedral_1_3412 = MDAnalysis.core.topologyobjects.Dihedral(
    GQ_1_3.select_atoms("name N1")+GQ_1_4.select_atoms("name N1")
    +GQ_1_1.select_atoms("name N1")+GQ_1_2.select_atoms("name N1")).value()

    dihedral_1_4123 = MDAnalysis.core.topologyobjects.Dihedral(
    GQ_1_4.select_atoms("name N1")+GQ_1_1.select_atoms("name N1")
    +GQ_1_2.select_atoms("name N1")+GQ_1_3.select_atoms("name N1")).value()
	
    

    dihedral_2_1234 = MDAnalysis.core.topologyobjects.Dihedral(
    GQ_2_1.select_atoms("name N1")+GQ_2_2.select_atoms("name N1")
    +GQ_2_3.select_atoms("name N1")+GQ_2_4.select_atoms("name N1")).value()

    dihedral_2_2341 = MDAnalysis.core.topologyobjects.Dihedral(
    GQ_2_2.select_atoms("name N1")+GQ_2_3.select_atoms("name N1")
    +GQ_2_4.select_atoms("name N1")+GQ_2_1.select_atoms("name N1")).value()

    dihedral_2_3412 = MDAnalysis.core.topologyobjects.Dihedral(
    GQ_2_3.select_atoms("name N1")+GQ_2_4.select_atoms("name N1")
    +GQ_2_1.select_atoms("name N1")+GQ_2_2.select_atoms("name N1")).value()

    dihedral_2_4123 = MDAnalysis.core.topologyobjects.Dihedral(
    GQ_2_4.select_atoms("name N1")+GQ_2_1.select_atoms("name N1")
    +GQ_2_2.select_atoms("name N1")+GQ_2_3.select_atoms("name N1")).value()
	
    dihedral[i] = [dihedral_1_1234,dihedral_1_2341,dihedral_1_3412,dihedral_1_4123,
    dihedral_2_1234,dihedral_2_2341,dihedral_2_3412,dihedral_2_4123]

    # The angles between the normals to the Gua planes and the axis of a quartet
    quartet_axis = ((GQ_1_1+GQ_1_2+GQ_1_3+GQ_1_4).center_of_mass()
    -(GQ_2_1+GQ_2_2+GQ_2_3+GQ_2_4).center_of_mass())/dist_GQ1_GQ2[i]
 
    r_N9N2 = (GQ_1_1.select_atoms("name N9").coordinates()
    -GQ_1_1.select_atoms("name N2").coordinates())[0]
    r_N9N2 = r_N9N2/np.linalg.norm(r_N9N2)
    r_N9O6 = (GQ_1_1.select_atoms("name N9").coordinates()
    -GQ_1_1.select_atoms("name O6").coordinates())[0]
    r_N9O6 = r_N9O6/np.linalg.norm(r_N9O6)
    norm_GQ_1_1 = np.cross(r_N9N2,r_N9O6)
    phi_quartet_axis_norm_GQ_1_1 = 180/np.pi*np.arccos(np.dot(quartet_axis,norm_GQ_1_1))

    r_N9N2 = (GQ_1_2.select_atoms("name N9").coordinates()
    -GQ_1_2.select_atoms("name N2").coordinates())[0]
    r_N9N2 = r_N9N2/np.linalg.norm(r_N9N2)
    r_N9O6 = (GQ_1_2.select_atoms("name N9").coordinates()
    -GQ_1_2.select_atoms("name O6").coordinates())[0]
    r_N9O6 = r_N9O6/np.linalg.norm(r_N9O6)
    norm_GQ_1_2 = np.cross(r_N9N2,r_N9O6)
    phi_quartet_axis_norm_GQ_1_2 = 180/np.pi*np.arccos(np.dot(quartet_axis,norm_GQ_1_2))

    r_N9N2 = (GQ_1_3.select_atoms("name N9").coordinates()
    -GQ_1_3.select_atoms("name N2").coordinates())[0]
    r_N9N2 = r_N9N2/np.linalg.norm(r_N9N2)
    r_N9O6 = (GQ_1_3.select_atoms("name N9").coordinates()
    -GQ_1_3.select_atoms("name O6").coordinates())[0]
    r_N9O6 = r_N9O6/np.linalg.norm(r_N9O6)
    norm_GQ_1_3 = np.cross(r_N9N2,r_N9O6)
    phi_quartet_axis_norm_GQ_1_3 = 180/np.pi*np.arccos(np.dot(quartet_axis,norm_GQ_1_3))

    r_N9N2 = (GQ_1_4.select_atoms("name N9").coordinates()
    -GQ_1_4.select_atoms("name N2").coordinates())[0]
    r_N9N2 = r_N9N2/np.linalg.norm(r_N9N2)
    r_N9O6 = (GQ_1_4.select_atoms("name N9").coordinates()
    -GQ_1_4.select_atoms("name O6").coordinates())[0]
    r_N9O6 = r_N9O6/np.linalg.norm(r_N9O6)
    norm_GQ_1_4 = np.cross(r_N9N2,r_N9O6)
    phi_quartet_axis_norm_GQ_1_4 = 180/np.pi*np.arccos(np.dot(quartet_axis,norm_GQ_1_4))

    r_N9N2 = (GQ_2_1.select_atoms("name N9").coordinates()
    -GQ_2_1.select_atoms("name N2").coordinates())[0]
    r_N9N2 = r_N9N2/np.linalg.norm(r_N9N2)
    r_N9O6 = (GQ_2_1.select_atoms("name N9").coordinates()
    -GQ_2_1.select_atoms("name O6").coordinates())[0]
    r_N9O6 = r_N9O6/np.linalg.norm(r_N9O6)
    norm_GQ_2_1 = np.cross(r_N9N2,r_N9O6)
    phi_quartet_axis_norm_GQ_2_1 = 180/np.pi*np.arccos(np.dot(quartet_axis,norm_GQ_2_1))

    r_N9N2 = (GQ_2_2.select_atoms("name N9").coordinates()
    -GQ_2_2.select_atoms("name N2").coordinates())[0]
    r_N9N2 = r_N9N2/np.linalg.norm(r_N9N2)
    r_N9O6 = (GQ_2_2.select_atoms("name N9").coordinates()
    -GQ_2_2.select_atoms("name O6").coordinates())[0]
    r_N9O6 = r_N9O6/np.linalg.norm(r_N9O6)
    norm_GQ_2_2 = np.cross(r_N9N2,r_N9O6)
    phi_quartet_axis_norm_GQ_2_2 = 180/np.pi*np.arccos(np.dot(quartet_axis,norm_GQ_2_2))

    r_N9N2 = (GQ_2_3.select_atoms("name N9").coordinates()
    -GQ_2_3.select_atoms("name N2").coordinates())[0]
    r_N9N2 = r_N9N2/np.linalg.norm(r_N9N2)
    r_N9O6 = (GQ_2_3.select_atoms("name N9").coordinates()
    -GQ_2_3.select_atoms("name O6").coordinates())[0]
    r_N9O6 = r_N9O6/np.linalg.norm(r_N9O6)
    norm_GQ_2_3 = np.cross(r_N9N2,r_N9O6)
    phi_quartet_axis_norm_GQ_2_3 = 180/np.pi*np.arccos(np.dot(quartet_axis,norm_GQ_2_3))

    r_N9N2 = (GQ_2_4.select_atoms("name N9").coordinates()
    -GQ_2_4.select_atoms("name N2").coordinates())[0]
    r_N9N2 = r_N9N2/np.linalg.norm(r_N9N2)
    r_N9O6 = (GQ_2_4.select_atoms("name N9").coordinates()
    -GQ_2_4.select_atoms("name O6").coordinates())[0]
    r_N9O6 = r_N9O6/np.linalg.norm(r_N9O6)
    norm_GQ_2_4 = np.cross(r_N9N2,r_N9O6)
    phi_quartet_axis_norm_GQ_2_4 = 180/np.pi*np.arccos(np.dot(quartet_axis,norm_GQ_2_4))

    phi_quartet_axis_norm_Gua[i] = [phi_quartet_axis_norm_GQ_1_1,
                                    phi_quartet_axis_norm_GQ_1_2,
                                    phi_quartet_axis_norm_GQ_1_3,
                                    phi_quartet_axis_norm_GQ_1_4,
                                    phi_quartet_axis_norm_GQ_2_1,
                                    phi_quartet_axis_norm_GQ_2_2,
                                    phi_quartet_axis_norm_GQ_2_3,
                                    phi_quartet_axis_norm_GQ_2_4]

    # The angles between Gua axes
    # the vector r_j_i joins the middle of the C4–C5 interval with the C8 atom in Gua_i in quartet j
    r_1_1 = (0.5*GQ_1_1.select_atoms("name C4").coordinates()+ 0.5*GQ_1_1.select_atoms("name C5").coordinates()
    -GQ_1_1.select_atoms("name C8").coordinates())[0]
    r_1_1 = r_1_1/np.linalg.norm(r_1_1)

    r_1_2 = (0.5*GQ_1_2.select_atoms("name C4").coordinates()+ 0.5*GQ_1_2.select_atoms("name C5").coordinates()
    -GQ_1_2.select_atoms("name C8").coordinates())[0]
    r_1_2 = r_1_2/np.linalg.norm(r_1_2)

    r_1_3 = (0.5*GQ_1_3.select_atoms("name C4").coordinates()+ 0.5*GQ_1_3.select_atoms("name C5").coordinates()
    -GQ_1_3.select_atoms("name C8").coordinates())[0]
    r_1_3 = r_1_3/np.linalg.norm(r_1_3)

    r_1_4 = (0.5*GQ_1_4.select_atoms("name C4").coordinates()+ 0.5*GQ_1_4.select_atoms("name C5").coordinates()
    -GQ_1_4.select_atoms("name C8").coordinates())[0]
    r_1_4 = r_1_4/np.linalg.norm(r_1_4)


    omega_1_12=180/np.pi*np.arccos(np.dot(r_1_1,r_1_2))
    omega_1_23=180/np.pi*np.arccos(np.dot(r_1_2,r_1_3))
    omega_1_34=180/np.pi*np.arccos(np.dot(r_1_3,r_1_4))
    omega_1_41=180/np.pi*np.arccos(np.dot(r_1_4,r_1_2))
    omega_1_13=180/np.pi*np.arccos(np.dot(r_1_1,r_1_3))
    omega_1_24=180/np.pi*np.arccos(np.dot(r_1_2,r_1_4))

    r_2_1 = (0.5*GQ_2_1.select_atoms("name C4").coordinates()+ 0.5*GQ_2_1.select_atoms("name C5").coordinates()
    -GQ_2_1.select_atoms("name C8").coordinates())[0]
    r_2_1 = r_2_1/np.linalg.norm(r_2_1)

    r_2_2 = (0.5*GQ_2_2.select_atoms("name C4").coordinates()+ 0.5*GQ_2_2.select_atoms("name C5").coordinates()
    -GQ_2_2.select_atoms("name C8").coordinates())[0]
    r_2_2 = r_2_2/np.linalg.norm(r_2_2)

    r_2_3 = (0.5*GQ_2_3.select_atoms("name C4").coordinates()+ 0.5*GQ_2_3.select_atoms("name C5").coordinates()
    -GQ_2_3.select_atoms("name C8").coordinates())[0]
    r_2_3 = r_2_3/np.linalg.norm(r_2_3)

    r_2_4 = (0.5*GQ_2_4.select_atoms("name C4").coordinates()+ 0.5*GQ_2_4.select_atoms("name C5").coordinates()
    -GQ_2_4.select_atoms("name C8").coordinates())[0]
    r_2_4 = r_2_4/np.linalg.norm(r_2_4)

    omega_2_12=180/np.pi*np.arccos(np.dot(r_2_1,r_2_2))
    omega_2_23=180/np.pi*np.arccos(np.dot(r_2_2,r_2_3))
    omega_2_34=180/np.pi*np.arccos(np.dot(r_2_3,r_2_4))
    omega_2_41=180/np.pi*np.arccos(np.dot(r_2_4,r_2_2))
    omega_2_13=180/np.pi*np.arccos(np.dot(r_2_1,r_2_3))
    omega_2_24=180/np.pi*np.arccos(np.dot(r_2_2,r_2_4))

    omega[i] = [omega_1_12,omega_1_23,omega_1_34,omega_1_41,omega_1_13,omega_1_24, 
                omega_2_12,omega_2_23,omega_2_34,omega_2_41,omega_2_13,omega_2_24]
    

np.savetxt('DistanceBetweenGuartets.dat',dist_GQ1_GQ2)
plt.figure()
plt.plot(dist_GQ1_GQ2) #plot only the second column
plt.xlabel('Frame')
plt.ylabel(r'Distance($\AA$)')
plt.title('Distance between the G-quartets')
#plt.style.use('bmh') #Uncomment to use matplotlib style sheets. Also, remove "import seaborn."
plt.savefig('DistanceBetweenGuartets.png',transparent=True)
plt.savefig('DistanceBetweenGuartets.pdf',transparent=True)

np.savetxt('DistanceBetweenComGuaToComGQuartets.dat',dist_GQ_Gua)
plt.figure()
plt.plot(dist_GQ_Gua) #plot only the second column
plt.xlabel('Frame')
plt.ylabel(r'Distance($\AA$)')
plt.title('Distance between COM of each Gua to the COM of G-quartets')
#plt.style.use('bmh') #Uncomment to use matplotlib style sheets. Also, remove "import seaborn."
plt.savefig('DistanceBetweenComGuaToComGQuartets.png',transparent=True)
plt.savefig('DistanceBetweenComGuaToComGQuartets.pdf',transparent=True)

np.savetxt('TwistAngle.dat',twistAngle)
plt.figure()
plt.plot(twistAngle) #plot only the second column
plt.xlabel('Frame')
plt.ylabel(r'Twist angle ($^{\circ}$)')
plt.title('Twist angle for quartets')
#plt.style.use('bmh') #Uncomment to use matplotlib style sheets. Also, remove "import seaborn."
plt.savefig('TwistAngle.png',transparent=True)
plt.savefig('TwistAngle.pdf',transparent=True)

np.savetxt('Dihedral.dat',dihedral)
plt.figure()
plt.plot(dihedral) #plot only the second column
plt.xlabel('Frame')
plt.ylabel(r'Dihedral angle ($^{\circ}$)')
plt.title('Dihedral angles')
#plt.style.use('bmh') #Uncomment to use matplotlib style sheets. Also, remove "import seaborn."
plt.savefig('Dihedral.png',transparent=True)
plt.savefig('Dihedral.pdf',transparent=True)

np.savetxt('Phi.dat',phi_quartet_axis_norm_Gua)
plt.figure()
plt.plot(phi_quartet_axis_norm_Gua) #plot only the second column
plt.xlabel('Frame')
plt.ylabel(r'Angles ($^{\circ}$)')
plt.title('Angles between the normals to the Gua planes and the axis of a quartet')
#plt.style.use('bmh') #Uncomment to use matplotlib style sheets. Also, remove "import seaborn."
plt.savefig('Phi.png',transparent=True)
plt.savefig('Phi.pdf',transparent=True)

np.savetxt('omega.dat',omega)
plt.figure()
plt.plot(omega) #plot only the second column
plt.xlabel('Frame')
plt.ylabel(r'Omega ($^{\circ}$)')
plt.title('Angles between Gua axes')
#plt.style.use('bmh') #Uncomment to use matplotlib style sheets. Also, remove "import seaborn."
plt.savefig('omega.png',transparent=True)
plt.savefig('omega.pdf',transparent=True)