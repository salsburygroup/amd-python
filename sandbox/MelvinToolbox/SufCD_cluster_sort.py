#!/usr/bin/env python

import numpy as np
import pandas as pd
import re
import os

dir = '/Volumes/RyanMdata/sufCandD/Docking/ClusterCluster/Output' #Where zdock outputs are

# Setup master data frame
models = pd.DataFrame(columns=["C","D","score"])

#Loop over all files in directory
for file in os.listdir(dir):
	if file.endswith(".out"):
		# Find out what c and d frame were used
		subunits = re.findall('DFrame([0-9]+)_m_CFrame([0-9]+)_m.',file)[0]
		DFrame = subunits[0]
		CFrame = subunits[1]
	
		# Import scores
		scores = np.loadtxt(dir+'/'+file,skiprows=5,usecols=(6,))
	
		# Make a temporary dataframe
		temp = pd.DataFrame({"C":"","D":"","score":scores})
		temp["C"]=CFrame
		temp["D"]=DFrame
		
		# Add temporary to master data frame
		models = pd.concat([models,temp])

models_sorted = models.sort(columns="score",ascending=False)

models.to_csv('/Volumes/RyanMdata/sufCandD/Docking/ClusterCluster/Output/models.csv',index=False)
models_sorted.to_csv('/Volumes/RyanMdata/sufCandD/Docking/ClusterCluster/Output/models_sorted.csv',index=False)
