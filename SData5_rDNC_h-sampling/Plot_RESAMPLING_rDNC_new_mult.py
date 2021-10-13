# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 11:38:54 2018
READS IN THE RESAMPLING OUTPUT AND CALCULATES STATISTICS FOR ALL RESAMPLING PARAMETERS
@author: Sasha
"""
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.gridspec as gridspec
from scipy import stats

spp = ['legitima', 'olangoensis', 'cingulifera']
#spp = ['pulex', 'laevis', 'melanica', 'longispina']
#spp = ['chaldaeus', 'ebraeus']#, 'sanguinolentus', 'miliaris']
#spp = ['tongmuensis', 'thomasi']

gene = 'COI'
maincolors = ['#E59866', '#BA4A00', '#6E2C00']#X/I
#maincolors = ['#A9DFBF', '#73C6B6', '#229954', '#145A32']#Daphnia
#maincolors = ['#1B4F72', '#2E86C1', '#85C1E9']#Conus
#maincolors = ['#A569BD', '#C39BD3'] #Tanytarsus

PVScolors = maincolors
TPVScolors = maincolors#['darkred', 'darkgreen', 'darkblue', 'darkgrey']
#TPVScolors = ['indianred', 'lime', 'dodgerblue', 'silver']

fig = plt.figure(figsize=(7, 4))#, constrained_layout=True)
spec = gridspec.GridSpec(ncols=1, nrows=1)#, figure=fig)
ax = fig.add_subplot(spec[0, 0])

ax.set_xlabel('Number of haplotypes sampled')
ax.set_ylabel('Priportion of EDV sDNCs')

for k in range(len(spp)):
    species = spp[k]
    #f = open('D:\sDNC_qHap_resampling\\'+species+'_COI.txt', 'r')
    f = open('D:\sDNC_qHap_refSPP_resampling\\'+species+'_'+gene+'.txt', 'r')
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    haps = []
    hap_values = {}
    for line in f:
        line.strip()
        data = line.split(' ')
        lineID = 'hap'+data[0]
        if int(data[0]) not in haps:
            haps.append(int(data[0]))
            hap_values[lineID] = [int(data[-1])]#number of species
        else:
            hap_values[lineID].append(int(data[-1]))
    #for i in range(len(spp_values['spp17'])):
    #     print len(spp_values['spp17'][i])
    f.close()
    #print '\n'+species
    hap_summ = []
    for n in haps:
        handle = 'hap'+str(n)
        hap_summ.append(float(sum(hap_values[handle]))/10)
        #print handle, hap_summ[-1]
    
    #Plots
    ax.errorbar(haps, hap_summ, linewidth=2, color=maincolors[k])
    
plt.legend(spp, loc="lower right")
plt.show()
#fig.savefig('E:\Python\Molecular_diagnoses\Resampling_output\Resampling_'+label+'-PVS.png')

