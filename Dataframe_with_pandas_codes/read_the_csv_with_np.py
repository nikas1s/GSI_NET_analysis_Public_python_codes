import numpy as np
import pandas as pd
from numpy import genfromtxt
data = genfromtxt('AME20_uncertainties.csv', delimiter=',')
#an example of reading a csv and printing the first 
print(data[:,1],data[:,2],data[:,3])
for i in range(0,len(data[:,1])):
    if int(data[:,1])==30 and data[:,2]==80:
        print (data[i][1],data[i][2],data[i][3])
