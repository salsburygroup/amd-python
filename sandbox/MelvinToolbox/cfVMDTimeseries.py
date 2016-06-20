import re
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import difflib

#Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description = 'Align to first frame and save new trajectory', add_help=False) 

inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-vmd', action='store', dest='vmd_cluster',help='VMD cluster output', type=str,required=True)
inputs.add_argument('-cf', action='store', dest='compare_cluster', help='Cluster timeseries to compare',required=True)

UserInput=parser.parse_args()

dfbig = pd.DataFrame(columns = ['frame', 'cluster_number'])

# Preprocess VMD input
with open (UserInput.vmd_cluster) as file:
    all_clusters = [result for result in re.findall('\{([\w+\s]+)',file.read(),re.S)] #re.S makes dot anything, even \n
with open (UserInput.vmd_cluster) as file:
    single_clusters = re.findall('\}\s([\w+\s]+)\s\{',file.read(),re.S)   

i = 1
for line in all_clusters:
    cluster = map(int, line.split(' '))
    dftemp = pd.DataFrame(columns = ['frame', 'cluster_number'])
    dftemp['frame'] = cluster
    dftemp['cluster_number'] = i
    dfbig = pd.concat([dfbig,dftemp])
    i = i+1

for line in single_clusters:
    cluster = map(int, line.split(' '))
    dftemp = pd.DataFrame(columns = ['frame', 'cluster_number'])
    dftemp['frame'] = cluster
    dftemp['cluster_number'] = i
    dfbig = pd.concat([dfbig,dftemp])
    i = i+1

dfsort = dfbig.sort_values(by=['frame'])
vmd_chain = dfsort['cluster_number'].values.astype(int)

# Match ClusterPrediction output scheme
vmd_chain = vmd_chain - 1
vmd_chain[np.where(vmd_chain == vmd_chain.max())] = -1

plt.figure()
plt.scatter(np.arange(len(vmd_chain)), vmd_chain, marker = '+')
plt.xlabel('Frame')
plt.ylabel('Cluster')
plt.title('VMD_QT')
plt.savefig(UserInput.vmd_cluster + '_timeseries.png')
plt.close()


# Get comparison input
cf_chain = np.genfromtxt(UserInput.compare_cluster)
cf_chain = cf_chain.astype(int)

sm = difflib.SequenceMatcher(None, vmd_chain, cf_chain)
print('Clustering identity is {}'.format(sm.ratio()))

