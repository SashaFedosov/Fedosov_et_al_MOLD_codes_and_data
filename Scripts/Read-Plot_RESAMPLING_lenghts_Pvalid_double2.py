# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 11:38:54 2018
READS IN THE RESAMPLING OUTPUT AND CALCULATES STATISTICS FOR ALL RESAMPLING PARAMETERS
@author: Sasha
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import stats
f = open("E:\Python\Molecular_diagnoses\Resampling_output\WALL_CLADES_resampling_Clade_Synallaxis_ND2-Barcode.txt", "r")
label = 'Synallaxis_ND2'
f.readline()
f.readline()
f.readline()
f.readline()
f.readline()
spp = []
spp_values = {}
spp_means = {}
spp_std = {}
spp_sem = {}
for line in f:
    line.strip()
    data = line.split(' ')
    if len(data) < 5:
        continue
    else:
        lineID = 'spp'+str(data[0])
        if lineID not in spp:
            spp.append(lineID)
            spp_values[lineID] = []
            spp_values[lineID].append([int(data[0])])#number of focus taxon species
            spp_values[lineID].append([int(data[1])])#number of focus taxon and non-focus taxa species
            spp_values[lineID].append([int(data[-1])])#combinations of 1 position found
        else:
            spp_values[lineID][0].append(int(data[0]))
            spp_values[lineID][1].append(int(data[1]))
            spp_values[lineID][2].append(int(data[-1]))#1 combinations found

Valid_means = []
Valid_std = []
Valid_sem = []
for n in spp:
    spp_means[n] = np.mean(spp_values[n][-1])
    Valid_means.append(np.mean(spp_values[n][-1]))
    spp_std[n] = np.std(spp_values[n][-1])
    Valid_std.append(np.std(spp_values[n][-1]))
    spp_sem[n] = stats.sem(spp_values[n][-1])
    Valid_sem.append(stats.sem(spp_values[n][-1]))

for o in spp:
    print o, spp_means[o]

g = open("E:\Python\Molecular_diagnoses\Resampling_output\WALL_CLADES_resampling_Clade_Synallaxis_ND2-Barcode2.txt", "r")
g.readline()
g.readline()
g.readline()
g.readline()
g.readline()
spp1 = []
spp_values1 = {}
spp_means1 = {}
spp_std1 = {}
spp_sem1 = {}
for line in g:
    line.strip()
    data = line.split(' ')
    if len(data) < 5:
        continue
    else:
        lineID = 'spp'+str(data[0])
        if lineID not in spp1:
            spp1.append(lineID)
            spp_values1[lineID] = []
            spp_values1[lineID].append([int(data[0])])#number of focus taxon species
            spp_values1[lineID].append([int(data[1])])#number of focus taxon and non-focus taxa species
            spp_values1[lineID].append([int(data[-1])])#combinations of 1 position found
        else:
            spp_values1[lineID][0].append(int(data[0]))
            spp_values1[lineID][1].append(int(data[1]))
            spp_values1[lineID][2].append(int(data[-1]))#1 combinations found

Valid_means1 = []
Valid_std1 = []
Valid_sem1 = []
for n in spp1:
    spp_means1[n] = np.mean(spp_values1[n][-1])
    Valid_means1.append(np.mean(spp_values1[n][-1]))
    spp_std1[n] = np.std(spp_values1[n][-1])
    Valid_std1.append(np.std(spp_values1[n][-1]))
    spp_sem1[n] = stats.sem(spp_values1[n][-1])
    Valid_sem1.append(stats.sem(spp_values1[n][-1]))

#Plots
fig, ax = plt.subplots(1, 1)
p1 = plt.errorbar([int(i[3:]) for i in spp], Valid_means, yerr = Valid_sem, color='xkcd:burnt orange')
p2 = plt.errorbar([int(i[3:]) for i in spp1], Valid_means1, yerr = Valid_sem1, color='xkcd:light orange')
#ax.grid()
#ax.set_xlabel('species sampled')
plt.xticks([int(i[3:]) for i in spp], rotation = 'vertical')
#ax.set_ylabel('Priportion of EDV-sDNCs')
ax.yaxis.set_ticks(np.arange(0, 1.2, 0.1))
#plt.title(label+' - EDV sDNCs')
ax.legend((p1, p2), ('Synallaxis', 'Synallaxis+4'), loc='lower right', shadow=True)
plt.show()
fig.savefig('E:\Python\Molecular_diagnoses\Resampling_output\Resampling_'+label+'-PVS.png')
"""
fig, ax = plt.subplots(1, 1)
plt.errorbar([int(i[3:]) for i in spp], Nind_mean, yerr = Nind_std, color='green')
ax.grid()
ax.set_xlabel('species sampled')
ax.set_ylabel('N independent diagnostic combinations')
plt.title(label+' - N Independent DNC')
fig.savefig('E:\Python\Molecular_diagnoses\Resampling_output\Resampling_'+label+'_2-Ind.png')

fig, ax = plt.subplots(1, 1)
p1 = plt.errorbar([int(i[3:]) for i in spp], One_found_mean, yerr = One_found_std, color='green')
p2 = plt.errorbar([int(i[3:]) for i in spp], Two_found_mean, yerr = Two_found_std, color='blue')
p3 = plt.errorbar([int(i[3:]) for i in spp], Three_found_mean, yerr = Three_found_std, color='orange')
p4 = plt.errorbar([int(i[3:]) for i in spp], Four_found_mean, yerr = Four_found_std, color='red')
ax.grid()
ax.legend((p1, p2, p3, p4), ('len(DNC) = 1', 'len(DNC) = 2', 'len(DNC) = 3', 'len(DNC) = 4'), loc='upper right', shadow=True)
ax.set_xlabel('species sampled')
ax.set_ylabel('N diagnostic combinations')
plt.title(label+' - All DNC by length')
plt.show()
fig.savefig('E:\Python\Molecular_diagnoses\Resampling_output\Resampling_'+label+'_3-All_D_combinations.png')

fig, ax = plt.subplots(1, 1)
plt.errorbar([int(i[3:]) for i in spp], One_valid_mean, yerr = One_valid_std, color='green')
plt.errorbar([int(i[3:]) for i in spp], Two_valid_mean, yerr = Two_valid_std, color='blue')
plt.errorbar([int(i[3:]) for i in spp], Three_valid_mean, yerr = Three_valid_std, color='orange')
plt.errorbar([int(i[3:]) for i in spp], Four_valid_mean, yerr = Four_valid_std, color='red')
ax.legend((p1, p2, p3, p4), ('len(DNC) = 1', 'len(DNC) = 2', 'len(DNC) = 3', 'len(DNC) = 4'), loc='upper left', shadow=True)
ax.grid()
ax.set_xlabel('species sampled')
ax.set_ylabel('N working diagnostic combinations')
plt.title(label+' - prop working DNC by length')
plt.show()
fig.savefig('E:\Python\Molecular_diagnoses\Resampling_output\Resampling_'+label+'_4-Proportion_Working_D_combinations.png')
"""
f.close()