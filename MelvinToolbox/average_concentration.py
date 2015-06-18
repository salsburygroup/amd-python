#! #!/usr/bin/env python

#Tue Jun  9 14:36:26 EDT 2015
#Ryan Melvin

#Script for getting a table of dssr value vs concentration for a given ion
#Unlike the matlab analog of this script, I've made the concentration list work in a sensible order

import re
from glob import glob
import pandas as pd

#Setup data frames
concentrations = [20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800]
columns= ['base pairs','multiplets','helices','stems','atom-base capping interactions','hairpin loops','bulges','internal loops', 'junctions', 'non-loop ss segments','kissing loops', 'A-minor motifs','ribose zippers','kink turns','phosphate interactions']
dssr_averages = pd.DataFrame(index=concentrations,columns=columns)

#Loop over each concentration and record average in correct place
for filename in glob("*/*/cat/dssr/dssr_values.txt"):
    concentration=re.findall('([0-9]+)mM',filename)
    concentration=map(int,concentration)
    df = pd.read_csv(filename)
    averages = df.mean()
    entry = averages.values
    dssr_averages.loc[concentration[0]]=entry

dssr_averages.to_csv('dssr_averages.txt')
dssr_averages.plot(subplots=True,layout=(5,3),figsize=(20,10),title='DSSR Averages')
plt.savefig('dssr_averages.png')
