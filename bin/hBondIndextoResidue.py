import argparse
import mdtraj as md
import pandas as pd

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description = 'Convert donor acceptor atom indecies and convert them to the appropriate residue', add_help=False) 

# List all possible user input
inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-pdb', action='store', dest='structure',help='Structure file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-hbonds', action='store', dest='hbonds',help='CSV of HBond occupancy data',type=str,required=True)
inputs.add_argument('-o', action='store', dest='outfile',help='File name for storing the result (for cytoscape use .txt)',type=str,required=True)

UserInput=parser.parse_args()

#import pdb file
structure = UserInput.structure
t=md.load(structure)
top=t.topology

#read in file with occupancy data
df = pd.read_csv(UserInput.hbonds, sep='\t', header=None)
donor = df[0]-1
acceptor = df[1]-1
occ = df[2]

i=0
#loop over first & second column and store residue
#print(donor[3])
f = open(UserInput.outfile, 'w')
while i < len(donor):
    
    donor_atom = top.atom(donor[i])
    acceptor_atom = top.atom(acceptor[i])
    occupancy = occ[i]
    
    print('%s%s %s%s %s' % (donor_atom.residue.name, donor_atom.residue.index+1, acceptor_atom.residue.name, acceptor_atom.residue.index+1, occupancy))
    f.write("{}{}\tpp\t{}{}\t{}\n".format(donor_atom.residue.name, donor_atom.residue.index+1, acceptor_atom.residue.name, acceptor_atom.residue.index+1, round(occupancy, 5)))
    i=i+1
    
