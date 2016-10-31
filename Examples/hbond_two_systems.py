#! /usr/bin/env python

# A custom hbond script made for comapring K and Na. Not sure if this is useful enough to make into a utility that
# obeys our API.
from __future__ import division
import MDAnalysis
import MDAnalysis.analysis.hbonds
import sys
import pandas as pd
import pickle
import numpy as np


K = MDAnalysis.Universe('/deac/salsburyGrp/melvrl13/Thrombin/FreeThrombinK/4dii_fill.B99990001_autopsf.pdb',
                        '/deac/salsburyGrp/melvrl13/Thrombin/FreeThrombinK/protein_all_aligned.dcd', permissive=True)
K_frames = len(K.trajectory)
K_dt = K.trajectory.dt

# Emulate VMD's detection algorithm if option is selected

    # Setup MD Analysis input
K_donors_acceptors = K.select_atoms(
    "name O or name O* or name O** or name N or name N* or name N** or name S or name S* or name S** or name F or name F* or name F**")
K_donors_acceptors_list = np.unique(K_donors_acceptors.names)


K_h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(
    K,
    'all',
    'all',
    selection1_type='both',
    update_selection1=False,
    update_selection2=False,
    distance=3.2,
    angle=120.0,
    donors=K_donors_acceptors_list,
    acceptors=K_donors_acceptors_list,
    detect_hydrogens='distance',
    distance_type='heavy'
)

K_h.run()
with open('K_h.pkl', 'w') as potassium:
    pickle.dump(K_h, potassium)

del K
K_h.generate_table()
K_df = pd.DataFrame.from_records(K_h.table)
del K_h
K_df['time']=K_df['time'].apply(lambda x: int(round(x / K_dt)))


Na = MDAnalysis.Universe('/deac/salsburyGrp/melvrl13/Thrombin/FreeThrombinNa/4dih_fill.B99990001_autopsf.pdb',
                         '/deac/salsburyGrp/melvrl13/Thrombin/FreeThrombinNa/protein_all_aligned.dcd', permissive=True)
Na_frames = len(Na.trajectory)
Na_dt = Na.trajectory.dt
Na_donors_acceptors = Na.select_atoms(
    "name O or name O* or name O** or name N or name N* or name N** or name S or name S* or name S** or name F or name F* or name F**")
Na_donors_acceptors_list = np.unique(Na_donors_acceptors.names)

Na_h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(
    Na,
    'all',
    'all',
    selection1_type='both',
    update_selection1=False,
    update_selection2=False,
    distance=3.2,
    angle=120.0,
    donors=Na_donors_acceptors_list,
    acceptors=Na_donors_acceptors_list,
    detect_hydrogens='distance',
    distance_type='heavy'
)

# Execute analysis.

Na_h.run()
with open('Na_h.pkl', 'w') as sodium:
    pickle.dump(Na_h, sodium)
del Na
Na_h.generate_table()
Na_df = pd.DataFrame.from_records(Na_h.table)
del Na_h
Na_df['time'] = Na_df['time'].apply(lambda x: 1 + max(K_df['time']) + int(round(x / Na_dt)))

both_df = pd.concat([K_df, Na_df])
del K_df
del Na_df

# Repeat for each frame individually and record
hbond_trajectory = pd.DataFrame(index=list(range(0, K_frames + Na_frames)))  # Empty data frame in memory
for frame in list(range(0, K_frames + Na_frames)):
    current_frame_hbonds = both_df.loc[
        both_df['time'] == int(frame), ['donor_resnm', 'donor_resid', 'acceptor_resnm', 'acceptor_resid']]
    current_frame_hbond_pairs = [row['donor_resnm'] + str(row['donor_resid']) + '-' + row['acceptor_resnm']
                                 + str(row['acceptor_resid']) for index, row in  current_frame_hbonds.iterrows()]
    for pair in current_frame_hbond_pairs:
        hbond_trajectory.loc[frame, pair] = 1  # Fill in the empties with 1 if the hbond occurs
    progress = "\r Motif calculation on Frame " + str(frame) + " of " + str(K_frames + Na_frames)  # status
    sys.stdout.write(progress)
    sys.stdout.flush()  # report status to terminal output
hbond_trajectory = hbond_trajectory.fillna(0)  # Fill any leftover empties with 0
hbond_trajectory.to_csv('trajectory.csv',index=False)  # Save as csv
