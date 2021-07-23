import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from math import ceil

#Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description = 'Histogram the clustering results', add_help=False) 

#List all possible user input
inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-title', action='store', dest='title',help='Title',type=str,required=True)
inputs.add_argument('-d', action='store', dest='data',help='Data',type=str,required=True)
inputs.add_argument('-m', action='store', dest='method',help='Clustering method',type=str,required=True)
inputs.add_argument('-sp', action='store', dest='sp',help='How many different simulations in the concatenated trajectory',type=int,required=True)
inputs.add_argument('-spnm', action='store', dest='spnm',help='Names for different simulations in the concatenated trajectory',type=str,required=True)

# Parse into useful form
UserInput = parser.parse_args()

data = pd.Series(np.loadtxt(UserInput.data))
delete_noise = "sed '/-1/'d " + "timeseries.txt" + " > " + "timeseries_no_noise.txt"
os.system(delete_noise)
data_no_noise=pd.Series(np.loadtxt('timeseries_no_noise.txt'))

## Graph one
fig1=plt.figure(1,figsize=(8,4))
data.plot.hist(bins=200, rwidth=0.9)
plt.title('Histogram of ' + UserInput.method + ' for ' + UserInput.title)
plt.xlabel('Clusters')
plt.ylabel('Number of frames')
fig1.savefig('Histogram_raw.png', dpi=200)
fig1.savefig('Histogram_raw.tiff', dpi=600)
fig1.savefig('Histogram_raw.pdf', bbox_inches='tight')

## Graph two
# Rearranging clusters
num_of_cluster=list(Counter(data_no_noise).values())
total_num_of_cluster_no_noise=len(data_no_noise)
total_num_of_cluster=len(data)
cluster_high2low=[]
sum_of_cluster_we_draw=[]

while np.max(num_of_cluster)/total_num_of_cluster_no_noise > 0.05:
    sum_of_cluster_we_draw.append(np.max(num_of_cluster))
    cluster_high2low.append(round(np.max(num_of_cluster)*100/total_num_of_cluster,2))
    num_of_cluster.remove(np.max(num_of_cluster))
    if len(num_of_cluster) == 0:
        break
rest=(len(data) - sum(sum_of_cluster_we_draw))
cluster_high2low.append(round(rest*100/total_num_of_cluster,2))

# X axis
index=[]
name_list = []
for i in range(len(cluster_high2low)-1):
    index = index + [i]
    name_list = name_list + [i]
index = index + [len(cluster_high2low)-1]
name_list = name_list + ['the rest']

# Y axis
Yindex=[]
uplimit=np.max(cluster_high2low)
if uplimit/10 > 5:
    Yname_list = [0]
    for i in range(ceil(uplimit/10)):
        Yindex = Yindex + [i*10]
        Yname_list = Yname_list + ['0.' + str(i+1)]
else:
    Yname_list = [0, 0.05]
    for i in range(ceil(uplimit/5)):
        Yindex = Yindex + [i*5]
        Yname_list = Yname_list + ['0.' + str((i+2)*5)]

# Plot
fig2=plt.figure(2,figsize=(8,4))
plt.title('Histogram of ' + UserInput.method + ' for ' + UserInput.title)
plt.xticks(index, name_list)
plt.yticks(Yindex, Yname_list)
plt.xlabel('Clusters')
plt.ylabel('Probability')
rects=plt.bar(range(len(cluster_high2low)), cluster_high2low)
for rect in rects:
    height = rect.get_height()
    plt.text(rect.get_x() + rect.get_width() / 2, height, str(height)+'%', ha='center', va='bottom')
fig2.savefig('histogram.png', dpi=200)
fig2.savefig('histogram.tiff', dpi=600)
fig2.savefig('histogram.pdf', bbox_inches='tight')

## Graph three
# Rearranging clusters
#interval=int(total_num_of_cluster/UserInput.sp)
WT=37131
dK9=26855
E8K=32850
R4A=37488

# interval1
interval1="sed -n '" + str(1) + "," + str(WT) + "p' timeseries.txt > timeseries1.txt"
os.system(interval1)

data1 = pd.Series(np.loadtxt('timeseries1.txt'))
if UserInput.method == 'HDBSCAN':
    delete_noise = "sed '/-1/'d " + "timeseries1.txt" + " > " + "timeseries_no_noise1.txt"
    os.system(delete_noise)
    data_no_noise1=pd.Series(np.loadtxt('timeseries_no_noise1.txt'))
else:
    data_no_noise1=pd.Series(np.loadtxt('timeseries1.txt')) 

# interval2
interval2="sed -n '" + str(WT+1) + "," + str(WT+dK9) + "p' timeseries.txt > timeseries2.txt"
os.system(interval2)
data2 = pd.Series(np.loadtxt('timeseries2.txt'))

if UserInput.method == 'HDBSCAN':
    delete_noise= "sed '/-1/'d " + 'timeseries2.txt' + " > " + 'timeseries_no_noise2.txt'
    os.system(delete_noise)
    data_no_noise2=pd.Series(np.loadtxt('timeseries_no_noise2.txt'))
else:
    data_no_noise2=pd.Series(np.loadtxt('timeseries2.txt')) 

# interval3
interval3="sed -n '" + str(WT+dK9+1) + "," + str(WT+dK9+E8K) + "p' timeseries.txt > timeseries3.txt"
os.system(interval3)

data3 = pd.Series(np.loadtxt('timeseries3.txt'))
if UserInput.method == 'HDBSCAN':
    delete_noise= "sed '/-1/'d " + 'timeseries3.txt' + " > " + 'timeseries_no_noise3.txt'
    os.system(delete_noise)
    data_no_noise3=pd.Series(np.loadtxt('timeseries_no_noise3.txt'))
else:
    data_no_noise3=pd.Series(np.loadtxt('timeseries3.txt')) 

# interval4
interval4="sed -n '" + str(WT+dK9+E8K+1) + "," + str(total_num_of_cluster) + "p' timeseries.txt > timeseries4.txt"
os.system(interval4)

data4 = pd.Series(np.loadtxt('timeseries4.txt'))
if UserInput.method == 'HDBSCAN':
    delete_noise= "sed '/-1/'d " + 'timeseries4.txt' + " > " + 'timeseries_no_noise4.txt'
    os.system(delete_noise)
    data_no_noise4=pd.Series(np.loadtxt('timeseries_no_noise4.txt'))
else:
    data_no_noise4=pd.Series(np.loadtxt('timeseries4.txt')) 

# Combine all simulations
if UserInput.method == 'HDBSCAN':
    separate_no_noise=[]
    for i in range(len(Counter(data))-1):    
        separate_no_noise.append([Counter(data1)[i], Counter(data2)[i], Counter(data3)[i], Counter(data4)[i]])
    separate=[]
    for i in range(len(Counter(data))):    
        separate.append([Counter(data1)[i-1], Counter(data2)[i-1], Counter(data3)[i-1], Counter(data4)[i-1]])
else:
    separate_no_noise=[]
    for i in range(len(Counter(data))):    
        separate_no_noise.append([Counter(data1)[i], Counter(data2)[i], Counter(data3)[i], Counter(data4)[i]])
    separate=separate_no_noise

# Pick clusters with big number of configurations
separate_no_noiseT=np.transpose(separate_no_noise)
#length=int(sum(list(Counter(data).values()))/UserInput.sp)
length=50000
separate_high2low=[]
sum_of_separate_high2low=[]

for i in range(len(separate_no_noiseT[0])):
    a=separate_no_noise[i]
    a_max=np.max(a)
    a_index=a.index(a_max)
    if a_max/sum(separate_no_noiseT[a_index]) > 0.1:
        sum_of_separate_high2low.append([a[0], a[1], a[2], a[3]])
        separate_high2low.append([round(a[0]*100/length,2), round(a[1]*100/length,2), round(a[2]*100/length,2), round(a[3]*100/length,2)])

# Order
if UserInput.method == 'HDBSCAN':
    def max2(elem):
        return np.max(elem)
    separate_high2low.sort(key=max2, reverse=True)
else:
    def WT(elem):
        return elem[0]
    separate_high2low.sort(key=WT, reverse=True)

# transpose
separate_high2lowT=np.transpose(separate_high2low)
sum_of_separate_high2lowT=np.transpose(sum_of_separate_high2low)

# Pick the rest and the noise
rest=([sum(separate_no_noiseT[0])-sum(sum_of_separate_high2lowT[0]), sum(separate_no_noiseT[1])-sum(sum_of_separate_high2lowT[1]), sum(separate_no_noiseT[2])-sum(sum_of_separate_high2lowT[2]), sum(separate_no_noiseT[3])-sum(sum_of_separate_high2lowT[3])])
separate_high2low.append([round(rest[0]*100/length,2), round(rest[1]*100/length,2), round(rest[2]*100/length,2), round(rest[3]*100/length,2)])

if UserInput.method == 'HDBSCAN':
    noise=(separate[0])
    separate_high2low.append([round(noise[0]*100/length,2), round(noise[1]*100/length,2), round(noise[2]*100/length,2), round(noise[3]*100/length,2)])

# transpose
separate_high2lowT2=np.transpose(separate_high2low)    

# X axis
bar_width=0.8/UserInput.sp
Xindex = np.arange(len(separate_high2lowT2[0]))
Xname_list = []
for i in range(len(separate_high2lowT[0])):
    Xname_list = Xname_list + [i]

Xname_list = Xname_list + ['rest']
if UserInput.method == 'HDBSCAN':
    Xname_list = Xname_list + ['noise']

# Y axis
Yindex=[]
uplimit=np.max(separate_high2lowT2)
if uplimit/10 > 5:
    Yname_list = [0]
    for i in range(ceil(uplimit/10)):
        Yindex = Yindex + [i*10]
        Yname_list = Yname_list + ['0.' + str(i+1)]
else:
    Yname_list = [0, 0.05]
    for i in range(ceil(uplimit/5)):
        Yindex = Yindex + [i*5]
        Yname_list = Yname_list + ['0.' + str((i+2)*5)]

# Plot
fig3=plt.figure(4,figsize=(16,4))
names=['WT', 'ΔK9', 'E8K', 'R4A']
rects1=plt.bar(Xindex + bar_width*0, separate_high2lowT2[0], bar_width, color='r', label=names[0])
rects2=plt.bar(Xindex + bar_width*1, separate_high2lowT2[1], bar_width, color='g', label=names[1])
rects3=plt.bar(Xindex + bar_width*2, separate_high2lowT2[2], bar_width, color='b', label=names[2])
rects4=plt.bar(Xindex + bar_width*3, separate_high2lowT2[3], bar_width, color='y', label=names[3])

plt.title('Histogram of ' + UserInput.method + ' for ' + UserInput.title)
plt.xticks(Xindex + bar_width*(UserInput.sp-1)/2, Xname_list)
plt.yticks(Yindex, Yname_list)
plt.xlabel('Clusters')
plt.ylabel('Probability')
plt.legend(loc='upper right')

def add_labels(rects):
    for rect in rects:
        height = rect.get_height()
        plt.text(rect.get_x() + rect.get_width() / 2, height, str(height)+'%', ha='center', va='bottom')
        rect.set_edgecolor('white')

add_labels(rects1)
add_labels(rects2)
add_labels(rects3)
add_labels(rects4)

fig3.savefig('separately_histogram.png', dpi=200)
fig3.savefig('separately_histogram.tiff', dpi=600)
fig3.savefig('separately_histogram.pdf', bbox_inches='tight')
