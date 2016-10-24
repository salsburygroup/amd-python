#! /usr/bin/env python

# A custom hbond script made for comapring F10 and FUMP10. Not sure if this is useful enough to make into a utility that
# obeys our API.
from __future__ import division
import MDAnalysis
import MDAnalysis.analysis.hbonds
import sys
import pandas as pd
import numpy as np


F10 = MDAnalysis.Universe('/Volumes/RyanMData/F10/F10.pdb', '/Volumes/RyanMData/F10/all5/CaKMgNaZn.dcd', permissive=True)
FUMP10 = MDAnalysis.Universe('/Volumes/RyanMData/FUMP10/FUMP10.pdb', '/Volumes/RyanMData/FUMP10/all5/CaKMgNaZn.dcd', permissive=True)
# Emulate VMD's detection algorithm if option is selected

    # Setup MD Analysis input
F10_donors_acceptors = F10.select_atoms(
    "name O or name O* or name O** or name N or name N* or name N** or name S or name S* or name S** or name F or name F* or name F**")
F10_donors_acceptors_list = np.unique(F10_donors_acceptors.names)

FUMP10_donors_acceptors = FUMP10.select_atoms(
    "name O or name O* or name O** or name N or name N* or name N** or name S or name S* or name S** or name F or name F* or name F**")
FUMP10_donors_acceptors_list = np.unique(FUMP10_donors_acceptors.names)

F10_h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(
    F10,
    'all',
    'all',
    selection1_type='both',
    update_selection1=False,
    update_selection2=False,
    distance=3.2,
    angle=120.0,
    donors=F10_donors_acceptors_list,
    acceptors=F10_donors_acceptors_list,
    detect_hydrogens='distance',
    distance_type='heavy'
)

FUMP10_h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(
    FUMP10,
    'not name O2',
    'not name O2',
    selection1_type='both',
    update_selection1=False,
    update_selection2=False,
    distance=3.2,
    angle=120.0,
    donors=FUMP10_donors_acceptors_list,
    acceptors=FUMP10_donors_acceptors_list,
    detect_hydrogens='distance',
    distance_type='heavy'
)

# Execute analysis.
F10_h.run()
FUMP10_h.run()
# Format as a table
F10_h.generate_table()
FUMP10_h.generate_table()
# Save as csv and change time to frames
F10_df = pd.DataFrame.from_records(F10_h.table)
F10_df['time']=F10_df['time'].apply(lambda x: int(round(x / F10.trajectory.dt))) #Convert from MDAnalysis default time values to frames
#F10_df.to_csv(UserInput.out_name + '_raw.csv', index=False) #Save as csv
FUMP10_df = pd.DataFrame.from_records(FUMP10_h.table)
FUMP10_df['time']=FUMP10_df['time'].apply(lambda x: max(F10_df['time']) + int(round(x / FUMP10.trajectory.dt)))

both_df = pd.concat([F10_df, FUMP10_df])

# Repeat for each frame individually and record
hbond_trajectory = pd.DataFrame(index=list(range(0, len(F10.trajectory) + len(FUMP10.trajectory)))) #Empty data frame in memory
for frame in list(range(0, len(F10.trajectory) + len(FUMP10.trajectory))):
    current_frame_hbonds = both_df.loc[both_df['time'] == int(frame), ['donor_resnm', 'donor_resid', 'acceptor_resnm', 'acceptor_resid']] #All hbonds in frame
    current_frame_hbond_pairs = [row['donor_resnm'] + str(row['donor_resid']) + '-' + row['acceptor_resnm'] + str(row['acceptor_resid'])
                                 for index, row in  current_frame_hbonds.iterrows()]
    for pair in current_frame_hbond_pairs:
        hbond_trajectory.loc[frame, pair] = 1 #Fill in the empties with 1 if the hbond occurs
    progress = "\r Motif calculation on Frame " + str(frame) + " of " + str(len(F10.trajectory) + len(FUMP10.trajectory)) #status
    sys.stdout.write(progress)
    sys.stdout.flush() #report status to terminal output
hbond_trajectory = hbond_trajectory.fillna(0) #Fill any leftover empties with 0
hbond_trajectory.to_csv('trajectory.csv',index=False) #Save as csv
