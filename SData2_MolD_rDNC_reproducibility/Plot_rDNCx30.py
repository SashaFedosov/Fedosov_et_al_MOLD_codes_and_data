# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 11:38:54 2018
THIS SCRIPT IN THE RESAMPLING OUTPUT (directory x30) AND PLOTS MEANS, SD AND min-max OF sDNC SCORING.
@author: Sasha
"""
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.gridspec as gridspec
from scipy import stats

spp = ['brasilensis','legitima', 'olangoensis', 'ebraeus']
#spp = ['pulex', 'laevis', 'melanica', 'longispina']
#spp = ['chaldaeus', 'ebraeus']#, 'sanguinolentus', 'miliaris']
#spp = ['tongmuensis', 'thomasi']

gene = 'COI'
maincolors = ['red', 'green', 'blue', 'black']
#maincolors = ['#E59866', '#BA4A00', '#6E2C00']#X/I
#maincolors = ['#A9DFBF', '#73C6B6', '#229954', '#145A32']#Daphnia
#maincolors = ['#1B4F72', '#2E86C1', '#85C1E9']#Conus
#maincolors = ['#A569BD', '#C39BD3'] #Tanytarsus

lists = {}
lists['legitima'] = ['[5:C]', '[5:C, 31:C]', '[5:C, 31:C, 70:C]', '[5:C, 31:C, 70:C, 133:T]']
lists['olangoensis'] = ['[403:A]', '[403:A, 421:A]', '[403:A, 421:A, 634:T]', '[403:A, 421:A, 634:T, 424:C]', '[403:A, 421:A, 634:T, 424:C, 418:C]']
lists['cingulifera'] = ['[265]', '[265, 466]', '[265, 466, 607]', '[265, 466, 607, 457]']
lists['ebraeus'] = ['[574:G, 586:A]', '[574:G, 586:A, 301:C]', '301:C, 550:T, 574:G, 586:A']
lists['brasilensis'] = ['[459:C]','[459:C, 648:A]','[459:C, 648:A, 653:T]','[459:C, 648:A, 653:T, 679:T]']

PVScolors = maincolors
TPVScolors = maincolors#['darkred', 'darkgreen', 'darkblue', 'darkgrey']
#TPVScolors = ['indianred', 'lime', 'dodgerblue', 'silver']

fig = plt.figure(figsize=(17, 4))#, constrained_layout=True)
spec = gridspec.GridSpec(ncols=4, nrows=1)#, figure=fig)
ax = fig.add_subplot(spec[0, 0])
ax1 = fig.add_subplot(spec[0, 1])
ax2 = fig.add_subplot(spec[0, 2])
ax3 = fig.add_subplot(spec[0, 3])
ax.set_xlabel('Number of positions in sDNC')
ax1.set_xlabel('Number of positions in sDNC')
ax2.set_xlabel('Number of positions in sDNC')
ax3.set_xlabel('Number of positions in sDNC')
ax.set_ylabel('sDNC score (out of 100)')

nums = [str(i) for i in range(10)[1:]]
for k in range(len(spp)):
    species = spp[k]
    #f = open('D:\sDNC_qHap_resampling\\'+species+'_COI.txt', 'r')
    f = open('D:\\'+species+'x30.sDNC.out', 'r')
    scores = {}
    for line in f:
        line.strip()
        data = line.split(' ')
        if data[0] in nums and int(data[0]) not in scores.keys():
            scores[int(data[0])] = [int(data[-1])]
        elif data[0] in nums and int(data[0]) in scores.keys():
            scores[int(data[0])].append(int(data[-1]))
        else:
            continue
    f.close()
    
    means = []
    std = []
    maxs = []
    mins = []
    print species
    for key in scores:
        means.append(np.mean(scores[key]))
        std.append(np.std(scores[key]))
        maxs.append(max(scores[key]))
        mins.append(min(scores[key]))
        print key, means[-1]
    #Plots
    #lens = [len(scores[key]) for key in scores]
    lens = lists[species]
    if k == 0:
        #ax.errorbar(scores.keys(), means, linewidth=3, yerr = std, elinewidth=1, alpha=0.5, color=maincolors[k])
        ax.xaxis.set_ticks(range(len(scores.keys())+1))
        ax.scatter(scores.keys(), means, color=maincolors[k])
        ax.errorbar(scores.keys(), means, linewidth=3, color=maincolors[k], yerr = std, elinewidth=1)
        for i, txt in enumerate(lens):
            ax.annotate(txt, (scores.keys()[i], means[i]+3))
        ax.errorbar(scores.keys(), maxs, linewidth=1, alpha=0.5, color=maincolors[k])
        ax.errorbar(scores.keys(), mins, linewidth=1, alpha=0.5, color=maincolors[k])
    elif k == 1:
        ax1.xaxis.set_ticks(range(len(scores.keys())+1))
        ax1.scatter(scores.keys(), means, color=maincolors[k])
        ax1.errorbar(scores.keys(), means, linewidth=3, color=maincolors[k], yerr = std, elinewidth=1)
        for i, txt in enumerate(lens):
            ax1.annotate(txt, (scores.keys()[i], means[i]+3))
        ax1.errorbar(scores.keys(), maxs, linewidth=1, alpha=0.5, color=maincolors[k])
        ax1.errorbar(scores.keys(), mins, linewidth=1, alpha=0.5, color=maincolors[k])
    elif k == 2:
        ax2.xaxis.set_ticks(range(len(scores.keys())+1))
        ax2.scatter(scores.keys(), means, color=maincolors[k])
        ax2.errorbar(scores.keys(), means, linewidth=3, color=maincolors[k], yerr = std, elinewidth=1)
        for i, txt in enumerate(lens):
            ax2.annotate(txt, (scores.keys()[i], means[i]+3))
        ax2.errorbar(scores.keys(), maxs, linewidth=1, alpha=0.5, color=maincolors[k])
        ax2.errorbar(scores.keys(), mins, linewidth=1, alpha=0.5, color=maincolors[k])
    elif k == 3:
        ax3.xaxis.set_ticks(range(len(scores.keys())+1))
        ax3.scatter(scores.keys(), means, color=maincolors[k])
        ax3.errorbar(scores.keys(), means, linewidth=3, color=maincolors[k], yerr = std, elinewidth=1)
        for i, txt in enumerate(lens):
            ax3.annotate(txt, (scores.keys()[i], means[i]+3))
        ax3.errorbar(scores.keys(), maxs, linewidth=1, alpha=0.5, color=maincolors[k])
        ax3.errorbar(scores.keys(), mins, linewidth=1, alpha=0.5, color=maincolors[k])
    else:
        print k

#plt.legend(spp, loc="lower right")
#plt.show()
#fig.savefig('E:\Python\Molecular_diagnoses\Resampling_output\Resampling_'+label+'-PVS.png')

