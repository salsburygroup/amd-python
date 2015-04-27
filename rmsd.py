#! /usr/bin/env python

#####
#rgyr.py
#Thu Apr 23 13:07:33 EDT 2015
#Ryan Melvin
#####
#Credit:http://www.mdanalysis.org/MDAnalysisTutorial/trajectories.html
#If you use this script, cite
#N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein. MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319-2327. doi:10.1002/jcc.21787
####
#TODO:
#####
#Example call
#python rgyr.py -structure '/Users/melvrl13/Documents/RyanM/F10/f10.psf' -traj '/Users/melvrl13/Desktop/foldingPaperWorking/weighted3200.dcd' -sel 'segid F10' -o '/Users/melvrl13/Desktop/testingPyplot/rgyr.dat' 

#Inputs:  structure file for trajectory (-structure) , trajectory (-traj), atomselection (-sel),  output name (-o)

#Outputs: rgyr time series in text file

#Dependencies
import MDAnalysis
import MDAnalysis.analysis.rms
import argparse
import numpy as np

#Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description = 'Calculate RGYR time series', add_help=False) 

#List all possible user input
inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-structure', action='store', dest='structure',help='Structure file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-ref', action='store', dest='ref_structure',help='Reference Structure',type=str,default=None)
inputs.add_argument('-refframe', action='store', dest='ref_frame',help='Reference Structure',type=int,default=0)
inputs.add_argument('-traj', action='store', dest='trajectory',help='Trajectory',type=str,required=True)
inputs.add_argument('-sel', action='store', dest='sel', help='Atom selection',type=str,default='all')
inputs.add_argument('-o', action='store', dest='out_name',help='Output file',type=str,required=True)

#Parse into useful form
UserInput=parser.parse_args()

#Define the universe (i.e., molecule in VMD)
u = MDAnalysis.Universe(UserInput.structure, 
        UserInput.trajectory
        )

if UserInput.ref_structure:
    reference_structure=MDAnalysis.Universe(UserInput.ref_structure)
else:
    reference_structure=None


#Set up calculation of RMSD
Rmsd = MDAnalysis.analysis.rms.RMSD(
        u, 
        reference=reference_structure, 
        select=UserInput.sel, 
        ref_frame=UserInput.ref_frame
        )

#Run calculation of RMSD
Rmsd.run()

#Save RMSD
save_rmsd=Rmsd.rmsd[:,2]
np.savetxt(UserInput.out_name,save_rmsd)
