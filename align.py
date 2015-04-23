#! /usr/bin/env python

#####
#AlignToPDB.py   
#Wed Apr 22 18:28:37 EDT 2015
#Ryan Melvin
#####
#Credit:
#If you use this script, cite
#Douglas L. Theobald (2005) "Rapid calculation of RMSD using a quaternion-based characteristic polynomial." Acta Crystallographica A 61(4):478-480.
#AND
#Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald (2010) "Fast determination of the optimal rotational matrix for macromolecular superpositions." J. Comput. Chem. 31, 1561-1563.
#AND
#N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein. MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319-2327. doi:10.1002/jcc.21787
####
#TODO:
#####



#Outputs: aligned dcd

#Example Call
#python align.py -structure '/Users/melvrl13/Documents/RyanM/F10/f10.psf' -traj '/Users/melvrl13/Desktop/foldingPaperWorking/weighted3200.dcd' -sel 'segid F10' -o '/Users/melvrl13/Desktop/testingPyplot/test.dcd'

import MDAnalysis
import MDAnalysis.analysis.align 
import argparse

#parse user input
#Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description = 'Align to first frame and save new trajectory', add_help=False) 

inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-structure', action='store', dest='structure',help='Structure file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-ref', action='store', dest='ref_structure',help='Reference Structure',type=str,default=None)
inputs.add_argument('-refframe', action='store', dest='ref_frame',help='Reference Structure',type=int,default=0)
inputs.add_argument('-traj', action='store', dest='trajectory',help='Trajectory',type=str,required=True)
inputs.add_argument('-sel', action='store', dest='sel', help='Atom selection',type=str,default='all')
inputs.add_argument('-o', action='store', dest='outName',help='Output file',type=str,required=True)

UserInput=parser.parse_args()

#Define the universe (i.e., molecule in VMD)
u = MDAnalysis.Universe(
        UserInput.structure, 
        UserInput.trajectory
        )

#Align according to user inputs
if UserInput.ref_structure:
    #Define the reference universe
    referencePDB = MDAnalysis.Universe(UserInput.ref_structure)
    MDAnalysis.analysis.align.rms_fit_trj(
            u, 
            referencePDB, 
            select=UserInput.sel, 
            filename=UserInput.outName
            )
else:
    #advance to reference frame
    u.trajectory[UserInput.ref_frame]
    MDAnalysis.analysis.align.rms_fit_trj(
            u, 
            u, 
            select=UserInput.sel, 
            filename=UserInput.outName
            )




