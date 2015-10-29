#! /usr/bin/env python

#####
#hbonds.py
#Thu Apr 23 15:43:40 EDT 2015
#Ryan Melvin
#####
#Credit:https://mdanalysis.googlecode.com/svn/trunk/doc/html/documentation_pages/analysis/hbonds.html
#If you use this script, cite
#N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein. MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319-2327. doi:10.1002/jcc.21787
#AND
#L.M. Gregoret, S.D. Rader, R.J. Fletterick, and F.E. Cohen. Hydrogen bonds involving sulfur atoms in proteins. Proteins, 9(2):99-107, 1991. 10.1002/prot.340090204.
####
#TODO:
#####
#Example call
#python /Users/melvrl13/Documents/AMD/AMD-PYTHON/Analysis/hbonds.py -structure /Volumes/RyanMdata/sufCandD/RelaxationSimulations/SufCD_with_ATP/Homology/LastFrameOfHomologyRound2/Complex.psf -t /Volumes/RyanMdata/sufCandD/RelaxationSimulations/SufCD_with_ATP/Homology/LastFrameOfHomologyRound2/round1/SufCD_ATP_relax_strip_stride.dcd -sel1 all -sel2 all -sel1_type both -d 4 -a 60 -o /Volumes/RyanMdata/sufCandD/RelaxationSimulations/SufCD_with_ATP/Homology/LastFrameOfHomologyRound2/round1/hbond_test.txt

#Outputs: Hydrogen bonds table with columns 0)time 1)donor index 2) acceptor index 3) donor residue name 4) donor residue id 5) donor atom 6) acceptor residue name 7) acceptor residue id 8) acceptor atom 9) distance 10) angle
#NOTE 1-based atom indexing
#NOTE 0-based column indexing

# Dependencies
from __future__ import division
import MDAnalysis
import MDAnalysis.analysis.hbonds
import argparse
import recfile
import fileinput
import sys
import pandas as pd


#Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
        description = (
            'output columns: 0)time 1)donor index 2) acceptor index 3) donor residue name 4) donor residue id 5) donor atom 6) acceptor residue name 7) acceptor residue id 8) acceptor atom 9) distance 10) angle'
            ), 
        add_help=False
        ) 

# List all possible user input
inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-structure', action='store', dest='structure',help='Structure file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-traj', action='store', dest='trajectory',help='Trajectory',type=str,required=True)
inputs.add_argument('-sel1', action='store', dest='sel1', help='Atom selection 1',type=str,default='protein')
inputs.add_argument('-sel2', action='store', dest='sel2', help='Atom selection 2',type=str,default='all')
inputs.add_argument('-sel1_type', action='store', dest='sel1_type', help='Selection 1 type, i.e. donor, acceptor or both',type=str,default='both')
inputs.add_argument('-update_sel1', action='store_true', dest='update_sel1', help='Upate selection 1 each frame?',default=False)
inputs.add_argument('-update_sel2', action='store_true', dest='update_sel2', help='Upate selection 2 each frame?',default=False)
inputs.add_argument('-d', '--distance', action='store', dest='distance',help='Distance cutoff in angstroms',type=float,default=3.0)
inputs.add_argument('-a', '--angle',  action='store', dest='angle', help='Angle cutoff in degrees', type=float, default=120.00)
inputs.add_argument('-o', action='store', dest='out_name',help='Output file',type=str,required=True)
inputs.add_argument('--extra_acceptors', action='store', dest='acceptors',help='Acceptors in addition to those defined by charmm27',type=str,default=None)
inputs.add_argument('--extra_donors', action='store', dest='donors',help='Donors in addition to those defined by charmm27',type=str,default=None)


# Parse into useful form
UserInput=parser.parse_args()

# Define the universe (i.e., molecule in VMD)
u = MDAnalysis.Universe(UserInput.structure, UserInput.trajectory, permissive=True)

# Setup analysis.
h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(
        u, 
        UserInput.sel1, 
        UserInput.sel2, 
        selection1_type=UserInput.sel1_type, 
        update_selection1=UserInput.update_sel1, 
        update_selection2=UserInput.update_sel2, 
        distance=UserInput.distance, 
        angle=UserInput.angle, 
        donors=UserInput.donors,
        acceptors=UserInput.acceptors,
        detect_hydrogens='distance'
        )

# Execute analysis.
h.run()

# Format as a table
h.generate_table()

# Save as csv
recf = recfile.Open(UserInput.out_name, 'w', delim=',')
header = 'time, donor_idx, acceptor_idx, donor_resnm, donor_resid, donor_atom, acceptor_resnm, acceptor_resid, acceptor_atom, distance, angle'
recf.fobj.write(header+'\n')
recf.write(h.table)
recf.close()

# Clean up the output file for humans and other csv-accepting software
for line in fileinput.input([UserInput.out_name],inplace=True):
    sys.stdout.write(line.replace('\00',''))

# Correct time from the MDAnalysis default (read stupid and infuriating) timestep
df = pd.read_csv(UserInput.out_name)
df['time']=df['time'].apply(lambda x: x/u.trajectory.dt)
df.to_csv(UserInput.out_name,index=False)
