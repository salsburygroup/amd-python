import pandas as pd
import os

path=os.getcwd()
n=10
m=6

whole=[]
for j in range(1,m):
    t=[]
    for i in range(n):
        t.append(pd.read_csv(path + '/' + str(j) + '/sub' + str(i) + '/hbond_trajectory.csv'))
    concat1 = pd.concat(t,ignore_index=True) 
    whole.append(concat1)
concat2 = pd.concat(t,ignore_index=True)
concat2 = concat2.fillna(0)
print(concat2)
concat2.to_csv('hbond_trajectory.csv', index=False)


#print(whole[0][0])
#print(whole[4][3])
#print(whole[3][9])
