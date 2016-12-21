import argparse
import mdtraj as md
import pandas as pd
import numpy as np

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description = 'Convert res1 res2 atom indecies and convert them to the appropriate residue', add_help=False) 

# List all possible user input
inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-pdb', action='store', dest='structure',help='Structure file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-sel', action='store', dest='sel',help='Atom selection',type=str, default='name CA')
inputs.add_argument('-corrs', action='store', dest='corrs',help='tab delimited file of correlation data',type=str,required=True)
inputs.add_argument('-o', action='store', dest='outfile',help='File name for storing the result (for cytoscape use .txt)',type=str,required=True)
inputs.add_argument('--resid', dest='resid', action='store_true', default='True')
inputs.add_argument('--no-resid', dest='resid', action='store_false')

UserInput=parser.parse_args()

#import pdb file
structure = UserInput.structure
t=md.load(structure)
top=t.topology
sel_Ind=top.select(UserInput.sel)+1
#print(sel_Ind)
flag=UserInput.resid

print(flag)

#read in file with occupancy data
df = pd.read_csv(UserInput.corrs, sep='\t', header=None)
res1 = df[0]-1
res2 = df[1]-1
occ = df[2]
#print(sel_Ind,res1)
#print(sel_Ind)
#print(res1)
i=0
j=0
#loop over first & second column and store residue


sel1 = pd.DataFrame()
sel2 = pd.DataFrame()

while i < len(sel_Ind):
    data1 = (res1.loc[res1==sel_Ind[i]])
    data2 = (res2.loc[res2==sel_Ind[i]])
    sel1=sel1.append(data1)
    sel2=sel2.append(data2)
    i=i+1


columns1=list(sel1.columns.values)
#print(columns1)
columns2=list(sel2.columns.values)

list3=[]
for e in columns2:
    if e in columns1:
        list3.append(e)
        
#print(list3)


f = open(UserInput.outfile, 'w')

for m in list3:
    sel1_atom = top.atom(res1[m])
    sel2_atom = top.atom(res2[m])
    occupancy = occ[m]
    if flag == True:
        print('%s%s %s%s %s' % (sel1_atom.residue.name, sel1_atom.residue.index+1, sel2_atom.residue.name, sel2_atom.residue.index+1, occupancy))
        f.write("{}{}\tpp\t{}{}\t{}\n".format(sel1_atom.residue.name, sel1_atom.residue.index+1, sel2_atom.residue.name, sel2_atom.residue.index+1, round(occupancy, 5)))
    else:
        print('%s %s %s' % (sel1_atom.index,  sel2_atom.index, occupancy))
        f.write("{}\t{}\t{}\n".format(sel1_atom.index, sel2_atom.index, round(occupancy, 5)))
    
    i=i+1
    