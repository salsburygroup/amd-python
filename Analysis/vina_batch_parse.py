#!/usr/bin/env Python

from __future__ import division
import glob
import linecache
import pandas as pd
import re
import seaborn #To make histogram look pretty
import numpy
import matplotlib.pyplot as plt

# Initialize data frame for storing models
df = pd.DataFrame(columns=['File','top_energy'])
i = 0 #Row number in dataframe

# Find all the pdbqt outputs
pattern = '*_out.pdbqt'

# Process those outputs
for file in glob.glob(pattern):
   line2 =  linecache.getline(file,2) #Get second line
   # Get energy from second line; trying to keep memory usage low
   df.loc[i,'top_energy'] = float(re.findall(r'REMARK.*(-\d+\.\d+)',line2)[0])

    # Store in data frame; identify by file name
   df.loc[i,'File'] = re.findall(r'(.*)_out.pdbqt',file)[0]
   i = i + 1

# Save data
df = df.sort(columns='top_energy',ascending=True)
df.to_csv('TopModels.csv', index=False)

# Save histogram of data
hist, bins = numpy.histogram(df.loc[:,'top_energy'], normed=True)

# Until I can find a better way...
width = (bins[1]-bins[0])
center = (bins[:-1]+bins[1:])/2
fig = plt.figure()
plt.bar(center, hist, align='center', width=width)
plt.xlabel('Free energy of binding')
plt.ylabel('Normalized probability')
plt.title('Vina Top Energy Outputs')
fig.savefig('TopEnergiesHistogram.png')
plt.clf()
