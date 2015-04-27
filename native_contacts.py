#! /usr/bin/env python

#####
#NativeContacts.py   
#Wed Apr 22 16:07:48 EDT 2015
#Ryan Melvin
#####
#Credit:
#If you use this script, cite
#N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein. MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319-2327. doi:10.1002/jcc.21787
#####
#TODO:
#####

#Inputs: reference structure (-ref), structure file for trajectory (-structure) , trajectory (-traj), atomselection (-sel), radius (-r) output name (-o)

#Outputs: time series, average contacts matrix (in a compressed file), time series plot, average contacts map plot

### Example call
#python NativeContacts.py -ref '/Volumes/RyanMdata/F10/Folding/weightedSims3200/hairpin.pdb' -structure '/Volumes/RyanMdata/F10/Folding/weightedSims3200/f10.psf' -traj '/Volumes/RyanMdata/F10/Folding/weightedSims3200/weighted3200.dcd' -sel 'name P*' -r 12 -o '/Users/melvrl13/Desktop/testingPyplot/contacts'

#NOTE: This function will overwrite already exiting files of the assigned name 

#Dependencies
import MDAnalysis
import MDAnalysis.analysis.contacts
import argparse
import matplotlib.pyplot

#parse user input
#Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description = 'Calculate and plot native contacts', add_help=False) 

inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-ref', action='store', dest='refStructure',help='Native Structure',type=str,required=True)
inputs.add_argument('-structure', action='store', dest='structure',help='Structure file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-traj', action='store', dest='trajectory',help='Trajectory',type=str,required=True)
inputs.add_argument('-sel', action='store', dest='sel', help='Atom selection',type=str,default='all')
inputs.add_argument('-r', action='store', dest='radius', help='Cut off distance for native contacts',type=float,required=True)
inputs.add_argument('-o', action='store', dest='outName',help='Output files prefix',type=str,required=True)

#Now put those values in a form other functions can easily read.
userInput=parser.parse_args()

#Format output file name
dataOut=userInput.outName+'.dat'
timeSeriesPlotOut=userInput.outName+'_TimeSeries.png'
averageContactsMapOut=userInput.outName+'_AverageContacts.png'

#Define the reference "universe." "Universe" is analogous to "molecule" in VMD
ref = MDAnalysis.Universe(userInput.refStructure)

#Select atoms for comparison. Note vmd-like selection language
refSelect=ref.selectAtoms(userInput.sel)

#Define the "universe" for comparison
u = MDAnalysis.Universe(userInput.structure, userInput.trajectory)

#Setup (not run) analysis
#ContactAnalysis1 compares one trajectory to one reference. 
CA1=MDAnalysis.analysis.contacts.ContactAnalysis1(
    u, 
    selection=userInput.sel, 
    refgroup=refSelect, 
    radius=userInput.radius, 
    outfile=dataOut)

#Run anslysis
#force=True overwrites files if they already exist. 
CA1.run(store=True,force=True)
    
#Plot time series of percent native contacts. #MDAnalysis calls this q_t
CA1.plot(filename=timeSeriesPlotOut)

#Plot matrix of average contacts
averageContactsMap=CA1.plot_qavg(filename=averageContactsMapOut)
