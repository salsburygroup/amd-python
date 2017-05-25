#! /usr/bin/env python
import argparse
import datetime
from os import listdir, getcwd
from os.path import isfile, join

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
    description='Generate bash script to align all DCD files in a directory to a common reference frame',
    add_help=False
)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument(
    '-sd',
    '--sourcedir',
    action='store',
    dest='sourceDir',
    help='Directory containing ref. and traj. files',
    type=str,
    default=getcwd()
)
inputs.add_argument(
    '-od',
    '--outputdir',
    action='store',
    dest='outputDir',
    help='Sub-directory for output trajectories',
    type=str,
    default=join(getcwd(), datetime.datetime.now().strftime("%Y-%m-%d_%H%M_align_batch"))
)
inputs.add_argument(
    '-os',
    '--outputsh',
    action='store',
    dest='outputScript',
    help='Generated script',
    type=str,
    default='align_batch.sh'
)
inputs.add_argument(
    '-e',
    '--env',
    action='store',
    dest='condaEnv',
    help='Conda environment to use',
    type=str
)
inputs.add_argument(
    '-p',
    '--postfix',
    action='store',
    dest='postfix',
    help='Postfix for aligned files',
    type=str,
    default='aligned'
)
inputs.add_argument(
    '-r',
    '--run',
    action='store',
    help='Run created batch file'
)
# Parameters to be passed to align.py
inputs.add_argument(
    '-sel',
    action='store',
    dest='sel',
    help='Atom selection',
    type=str,
    default='all'
)


# Parse into useful form
UserInput = parser.parse_args()

#check for traj pdb then psf then error
onlyfiles = [f for f in listdir(UserInput.sourceDir) if isfile(join(UserInput.sourceDir, f))]

refs = [file for file in onlyfiles if file.endswith(".pdb") or file.endswith(".psf")]
trajs = [file for file in onlyfiles if file.endswith(".dcd")]

if len(refs) != 1:
    print("ERROR: Must be exactly 1 pdb or 1 psf present. Exiting!", file=sys.stderr)
    sys.exit()

if len(trajs) == 0:
    print("ERROR: No dcd found. Exiting!", file=sys.stderr)
    sys.exit()

#for each dcd, fold comprehension with +\n
def alignCmd(ref, traj, sel, outtraj):
    return """python align.py -ref "{0}" -traj "{1}" -sel "{2}" -o "{3}\"""".format(ref, traj, sel, outtraj)

body = ["#!/bin/bash"]
if not (UserInput.condaEnv is None):
    body.append("source activate " + UserInput.condaEnv)
else:
    print("WARNING: No Conda environment specified.")
body.append("mkdir \"" + UserInput.outputDir + "\"")

for traj in sorted(trajs):
    body.append(alignCmd(join(UserInput.sourceDir, refs[0]), join(UserInput.sourceDir, traj), UserInput.sel, join(UserInput.outputDir, traj[:-4] + "_" + UserInput.postfix + ".dcd")))

outfilePath = join(UserInput.sourceDir, UserInput.outputScript)
outfile = open(outfilePath, 'w')
outfile.write("\n".join(body))
outfile.close()

print("SUCCESS: Created \"" + outfilePath + "\". Exiting normally.")
#write body + dcd output

#run script from current directory if toggled

