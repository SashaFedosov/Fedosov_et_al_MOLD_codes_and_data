# -*- coding: utf-8 -*-
"""
This script produces plots from the pDNC resampling results.
Plots built for all pDNCs regardless of length, and quotes need to be adjusted to switch from robustness curve to scatterplot.
Status: This is the final script and does not require revisions.

"""
import sys, os
import numpy as np
from numpy.polynomial.polynomial import polyfit
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import stats
from scipy.stats import ttest_ind
sp = 'cingulifera'
checked_length = 3
f = open('D:\qHaps_pDNC_resampling\\'+sp+'.txt', 'r')

#label = 'laevis_COI'

f.readline()
f.readline()
f.readline()
f.readline()
hap_values = {}
datalist = []
for line in f:
    line = line.rstrip().replace(',','')
    things = line.split('|')
    data = [thing.split(' ') for thing in things]
    
    allone = sum([int(j) for j in data[1]])#4 proportion of pure EDVpDNCs
    if allone == 0:
        oneEDV = 0
    else:
        oneEDV = (float(data[1][0])+float(data[1][2]))/allone
                 
    alltwo = sum([int(j) for j in data[2]])#4 proportion of pure EDVpDNCs
    twoEDV = (float(data[2][0])+float(data[2][2]))/alltwo

    allthree = sum([int(j) for j in data[3]])#4 proportion of pure EDVpDNCs
    threeEDV = (float(data[3][0])+float(data[3][2]))/allthree

    allfour = sum([int(j) for j in data[4]])#4 proportion of pure EDVpDNCs
    fourEDV = (float(data[4][0])+float(data[4][2]))/allfour
    
    datalist.append([int(data[0][0]), float(data[0][3]), oneEDV, twoEDV, threeEDV, fourEDV])
f.close()

#Scatterplot fig3E
#title = sp+' - '+str(checked_length)+ ' pDNC diagnosis robustness - qPVS'

qPVS = [thing[1] for thing in datalist]
EDV1 = [thing[2] for thing in datalist]
EDV2 = [thing[3] for thing in datalist]
EDV3 = [thing[4] for thing in datalist]
EDV4 = [thing[5] for thing in datalist]

plt.scatter(qPVS, EDV1, color='blue', alpha=0.25)
plt.scatter(qPVS, EDV2, color='green', alpha=0.25)
plt.scatter(qPVS, EDV3, color='orange', alpha=0.25)

plt.yticks(np.arange(0, 1.2, 0.2))
plt.title (sp)
plt.xlabel('qPVS')
plt.ylabel('proportion of EDV-pDNCs')
plt.legend()
plt.show()
#plt.savefig('D:\\'+title+'.svg')

