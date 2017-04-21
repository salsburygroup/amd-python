#! /usr/bin/env python

# An example script to combine h-bond trajectory via pandas with sparse representation of data


import pandas as pd
#import pickle



df_1 = pd.read_csv('/home/xiaoj12/freeThrombin/K/run1/hbond_trajectory.csv')
df_2 = pd.read_csv('/home/xiaoj12/freeThrombin/K/run2/hbond_trajectory.csv')
df_3 = pd.read_csv('/home/xiaoj12/freeThrombin/K/run3/hbond_trajectory.csv')
df_4 = pd.read_csv('/home/xiaoj12/freeThrombin/K/run4/hbond_trajectory.csv')
df_5 = pd.read_csv('/home/xiaoj12/freeThrombin/K/run5/hbond_trajectory.csv')
K_df = pd.concat([df_1, df_2, df_3, df_4, df_5],ignore_index=True)

df_1 = pd.read_csv('/home/xiaoj12/freeThrombin/Na/run1/hbond_trajectory.csv')
df_2 = pd.read_csv('/home/xiaoj12/freeThrombin/Na/run2/hbond_trajectory.csv')
df_3 = pd.read_csv('/home/xiaoj12/freeThrombin/Na/run3/hbond_trajectory.csv')
df_4 = pd.read_csv('/home/xiaoj12/freeThrombin/Na/run4/hbond_trajectory.csv')
df_5 = pd.read_csv('/home/xiaoj12/freeThrombin/Na/run5/hbond_trajectory.csv')
Na_df = pd.concat([df_1, df_2, df_3, df_4, df_5],ignore_index=True)

del df_1, df_2, df_3, df_4, df_5
both_df = pd.concat([K_df, Na_df], ignore_index=True)


both_df = both_df.fillna(0)  # Fill any leftover empties with 0
both_df.to_csv('trajectory_KNa.csv',index=False)  # Save as csv
