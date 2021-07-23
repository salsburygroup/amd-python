import numpy as np
import os

D=[{1,2},{3,4},{5,6,}]
F=[]

for i in range(len(D)):
    F.append(list(D[i]))
#np.savetxt('/home/wud18/python/wdz.txt', F, fmt='%d', delimiter=',')
for i in range(len(D)):
    np.savetxt('/home/wud18/python/wdz'+str(i)+'.txt', list(D[i]), fmt='%i')
#a=np.loadtxt('/home/wud18/python/wdz.txt', fmt='%s')
#print(a)

#a=np.loadtxt('/home/wud18/python/wdz.txt')



#for i in range(len(D)):


