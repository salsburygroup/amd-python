#!/usr/bin/env Python

import pandas as pd
import csv
import re
import os

# Initialize data frame for storing models
df = pd.DataFrame(columns=['cluster','cluster_population','model','energy','in_site'])
i = 0 #Row number in dataframe

# Get cluster populations
with open('/Volumes/RyanMdata/F10/Albumin/F10ClusterPopulations.txt','rb') as populations:
    reader=csv.reader(populations, quoting=csv.QUOTE_NONNUMERIC)
    rows=list(reader)
    pop = [float(j[0]) for j in rows]
for cluster in range(1,52):
    dir='/Volumes/RyanMdata/F10/Albumin/ensemble_longAcid/output/cluster%s' %cluster
    acid_sites = dir + '/acid_sites.txt'
    population=pop[cluster - 1]
    for file in os.listdir(dir):
        if re.match("Frame\d+_out.txt",file):
            model_details = dir+ '/' + file
    # Get list of models with atoms in acid site
    with open(acid_sites,'rb') as f:
        reader = csv.reader(f)
        sites = list(reader)[0] #record the entries in text file as list
        # Get model # and energy for each model
    with open(model_details,'rb') as f2:
        for line in f2: #Loop over lines in each file
            model = re.findall(r'\s+([1-9])\s+(-[0-9]+\.[0-9]+)',line) #creates a list of single tuple [(model,engergy0]
            if model:
                model = list(model[0]) #convert tuple to list
                # Put cluster number and population at beginning of list
                model.insert(0,cluster)
                model.insert(1,population)
                # Check if this model is in the acid site list
                if model[2] in sites:
                    model.append(1) #True
                else:
                    model.append(0) #False
                #Add line in dataframe for model
                df.loc[i]=model 
                i=i+1



