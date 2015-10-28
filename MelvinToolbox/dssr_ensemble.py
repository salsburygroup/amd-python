import argparse
from MDAnalysis import Universe
from split_frames import split_all
from os import path, makedirs, system
from re import findall
from pandas import DataFrame


parser = argparse.ArgumentParser(
        description = (
            'run dssr on a trajectory'
            ), 
        add_help=False
        ) 

inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-s', action='store', dest='structure',help='Structure file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-t', action='store', dest='trajectory',help='Trajectory',type=str,required=True)
inputs.add_argument('-o', action='store', dest='out_dir',help='Output file',type=str,required=True)

UserInput=parser.parse_args()

working_directory = UserInput.out_dir
structure = UserInput.structure
trajectory = UserInput.trajectory

#make a directory for temporarily storing frame PDBs
frame_directory = working_directory + "/Frames"
if not path.exists(frame_directory):
    makedirs(frame_directory)


#Get number of frames
u = Universe(structure, trajectory)
frames = len(u.trajectory)
u = None

#Split DCD into frame PDBs
frame_name = frame_directory + "/Frames"
split_all(structure,trajectory,frame_name)

#Make directories for storing dssr output
terminal_output_directory = working_directory + "/terminal_out"
po4_output_directory = working_directory + "/po4_out"
if not path.exists(terminal_output_directory):
    makedirs(terminal_output_directory)
if not path.exists(po4_output_directory):
    makedirs(po4_output_directory)

#Run DSSR on each pdb and record the po4 and terminal outputs to text files
for j in range(frames):
    terminal_output_name=terminal_output_directory + "/Frame" + str(j) + ".txt"
    command='/Applications/x3DNA/x3dna-dssr --po4 -i=' + frame_name + str(j) + '.pdb -o=' + po4_output_directory + '/Frame' + str(j) + '.txt >&' + terminal_output_name         
    #args=split(command)
    #with open(terminal_output_name,"wb") as log:
        #process = Popen(args, stdout=log, cwd=working_directory)
    system(command) 
    #This is sloppy. An attempt at the "right" way is commented above, but output consistently goes to the terminal window rather than the file when I use it

#Get the values from each iteration and record them in a pandas data frame
#Setup dataframe
columns= ['base pairs','multiplets','helices','stems','atom-base capping interactions','hairpin loops','bulges','internal loops', 'junctions', 'non-loop ss segments','kissing loops', 'A-minor motifs','ribose zippers','kink turns','phosphate interactions']
df = DataFrame(columns=columns)
df.index.name='Frame'

#Get values and write to data frame
for j in range(frames):
    terminal_output_name=terminal_output_directory + "/Frame" + str(j) + ".txt"
    report = open(terminal_output_name)
    values = []
    for line in report:
        line=line.rstrip()
        x=None
        x=findall('^total.*: ([0-9]+)',line)
        if len(x)>0:
            values.append(x[0])
    values = map(int,values)
    df.loc[j] = values

#Save data frame as csv for later plotting or manipulation
csv_name = working_directory + '/dssr_values.txt'
df.to_csv(csv_name,index=False)
            
