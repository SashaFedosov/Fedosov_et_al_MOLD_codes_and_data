#This script finds a characteristic combination of nucleotide positions for each clade
import os, sys
import random
#1. FUNCTIONS OF THE MolD##################################################################################################################################################
#1.1. SORTING ENTRIES BY CLADE AND IDENTIFYING NUCLEOTIDE POSITIONS SHARED WITHIN CLADE
def Step1(raw_records):
    Clades=[]
    for i in range(len(raw_records)):
        Clade=raw_records[i][1]
        if Clade not in Clades:
            Clades.append(Clade)
    clade_sorted_seqs = {}
    clade_sorted_spp = {}
    for letter in Clades:
        clade_sorted_seqs[letter]=[]
        clade_sorted_spp[letter]=[]
        for i in range(len(raw_records)):
            if raw_records[i][1]==letter:
                clade_sorted_seqs[letter].append(raw_records[i][2])
                if raw_records[i][0] not in clade_sorted_spp[letter]:
                    clade_sorted_spp[letter].append(raw_records[i][0])
    shared_positions={}
    for key in clade_sorted_seqs:
        sh_pos=[]
        for i in range(len(clade_sorted_seqs[key][0])):
            shared_nucleotide = True
            csm = clade_sorted_seqs[key][0][i] #candidate shared nucleotide
            for j in range(1, len(clade_sorted_seqs[key])):
                if clade_sorted_seqs[key][j][i] != csm:
                    shared_nucleotide = False
                    break
            if shared_nucleotide == True:
                sh_pos.append(i)
        shared_positions[key]=sh_pos
    return Clades, clade_sorted_seqs, shared_positions, clade_sorted_spp

#1.2. COMPILING COMPARISON LISTS FOR CLADES AND IDENTIFYING VARIABLE POSITIONS AND N PRIORITY POSITIONS WITH LARGEST CUTOFFS
def C_VP_PP(clade_sorted_seqs, clade, shared_positions, CUTOFF):# complist_variable_positions_priority_positions; Arguments: dictionary, string, dictionary
    CShN={}
    for pos in shared_positions[clade]:
        CShN[pos] = clade_sorted_seqs[clade][0][pos]#creates a dictionary shared position : nucleotide
    complist=[]
    for key in clade_sorted_seqs:
        if key != clade:
            complist = complist + clade_sorted_seqs[key]#creates a list of all other sequences for comparison
#    cstrings={}
#    for key in CShN:
#        cstrings[key] = ''
#        for j in complist:
#            cstrings[key] = cstrings[key] + j[key]
#    variable_positions = {}
#    for key in cstrings:
#        if cstrings[key].count(CShN[key]) < len(cstrings[key]):
#            variable_positions[key] = cstrings[key]
    cutoffs = {}
    for key in CShN:
        newcomplist = []
        for k in complist:
            if k[key] == CShN[key]:
                newcomplist.append(k)
            else: continue
        cutoffs[key] = len(complist) - len(newcomplist)
    CPP = []
    for key in sorted(cutoffs, key = cutoffs.get, reverse = True):
        CPP.append(key)
    Clade_priority_positions = {}
    for position in CPP[:CUTOFF]:#Here you define how many of the clade shared combinations are used in subsequent search
        Clade_priority_positions[position] = CShN[position]
    return complist, Clade_priority_positions, CShN

#1.3 RANDOM SEARCH ACROSS PRIORITY POSITIONS TO FIND AND REFINE CANDIDATE DIAGNOSTIC COMBINATIONS
def random_position(somelist, checklist):#gives a random index (integer) of the specified range, and returns indexed somelist element if it is not present in the checklist 
    while True:
        i = random.randint(0, len(somelist) - 1)
        if somelist[i] not in checklist:
            return somelist[i]
            break
        else:
            continue

def step_reduction_complist(clade, complist, CPP, checked_ind):#checks randomly selected positions of CladeSharedNucleotides with sequences of other clades, until a diagnostic combination of nucleotides for a selected clade is found.
    if len(complist) == 0:
        return checked_ind
    elif len(checked_ind) == len(CPP):
        return checked_ind
    else:
        newcomplist = []
        pos = random_position(list(CPP.keys()), checked_ind)
        for j in complist:
            if j[pos] == CPP[pos]:
                newcomplist.append(j)
            else: continue
        new_checked_ind = checked_ind + [pos]
        return step_reduction_complist(clade, newcomplist, CPP, new_checked_ind)
def ConditionD(newcomb, complist, CPP):#The function checks the 'Condition D' - i.e. whither any given combination of nucleotide positions is diagnostic for the selected clade
    ContD = False
    for i in newcomb:
        newcomplist = []
        for m in complist:
            if m[i] == CPP[i]:
                newcomplist.append(m)
            else: continue
        complist = newcomplist
    if len(complist) == 0:
        ContD = True
    return ContD
def RemoveRedundantPositions(raw_comb, complist, CPP):# The function removes positions from the raw combinations one by one, and then checks whether new combination fulfills the condition D, thus recursively reducing the diagnostic combination.
    red_possible = False
    for j in raw_comb:
        newcomb = [k for k in raw_comb if k != j]
        if ConditionD(newcomb, complist, CPP) == True:
            red_possible = True
            return RemoveRedundantPositions(newcomb, complist, CPP)
        else: pass
    if red_possible == False:
        return raw_comb

#1.4 IMPLEMENTS PREVIOUS FUNCTIONS TO FIRST PERFORM THE ITERATED RANDOM SEARCH, THEN REFINE OBTAINED COMBINATIONS, AND THEN RETURN MOLECULAR DIAGNOSIS
def Diagnostic_combinations(qCLADE, complist, CPP, n1, maxlen1):
    Achecked_ind = []
    bestlists = []
    n = n1
    while n>0:#STEP3 proposes raw diagnostic combinations
        m = step_reduction_complist(qCLADE, complist, CPP, Achecked_ind)
        if len(m) < maxlen1 and sorted(m) not in bestlists:
            bestlists.append(sorted(m))
        n=n-1
    bestlists.sort(key=len)
    priority_lists = bestlists[:500]
    #print len(priority_lists), len(priority_lists[0]), len(priority_lists[-1])
    diagnostic_combinations = []
    for raw_comb in priority_lists:
        refined_comb = RemoveRedundantPositions(raw_comb, complist, CPP)
        if not refined_comb in diagnostic_combinations:
            diagnostic_combinations.append(refined_comb)
    diagnostic_combinations.sort(key=len)
    return diagnostic_combinations

def IndependentKey(diagnostic_combinations):#
    independent_combinations = []
    selected_positions = []
    for i in range(len(diagnostic_combinations)):
        if len(selected_positions) == 0:
            for j in range(0, i):
                if len(set(diagnostic_combinations[i]) & set(diagnostic_combinations[j])) == 0 and len(set(diagnostic_combinations[i]) & set(selected_positions)) == 0:
                    independent_combinations.append(diagnostic_combinations[i])
                    independent_combinations.append(diagnostic_combinations[j])
                    for k in range(len(diagnostic_combinations[i])):
                        selected_positions.append(diagnostic_combinations[i][k])
                    for l in range(len(diagnostic_combinations[j])):
                        selected_positions.append(diagnostic_combinations[j][l])
        else:
            if len(set(diagnostic_combinations[i]) & set(selected_positions)) == 0:
                independent_combinations.append(diagnostic_combinations[i])
                for k in range(len(diagnostic_combinations[i])):
                    selected_positions.append(diagnostic_combinations[i][k])
    key_positions = []
    for pos in diagnostic_combinations[0]:
        KP = True
        for combination in diagnostic_combinations[1:]:
            if pos not in combination:
                KP = False
                break
            else: continue
        if KP == True:
            key_positions.append(pos)
    return independent_combinations, key_positions

#1.5 FUNCTIONS TO BUILD AND REFINE sDNC
def random_sequence(SEQ, Pdiff):#SEQ is a string and the Pdiff is an integer
    N = int(len(SEQ)*Pdiff/100)
    PosToChange = random.sample(range(0, len(SEQ)), N)
    Nucleotides = ['A', 'C', 'G', 'T']
    NEWSEQ = ''
    for i in range(len(SEQ)):
        if i not in PosToChange:
            NEWSEQ = NEWSEQ + SEQ[i]
        else:
            NEWSEQ = NEWSEQ + random.sample([j for j in Nucleotides if j != SEQ[i]], 1)[0]
    return NEWSEQ

def GenerateBarcode(Diagnostic_combinations, length):#This function calculates diagnostic combinations and assembles a barcode of desired length for a query taxon
    positions_involved = []
    positions_freq={}
    for comb in Diagnostic_combinations:
        if len(comb) > 1:
            for i in comb:
                positions_involved.append(i)
    for pos in set(positions_involved):
        positions_freq[pos] = positions_involved.count(pos)#This part needs to be re-written to first include the single-nucleotide combinations.
    barcode = [i[0] for i in Diagnostic_combinations if len(i) == 1] + [i for i in sorted(positions_freq, key=positions_freq.get, reverse=True)]
    return barcode[:length]

def GenerateBarcode2(Diagnostic_combinations, length):#This function calculates diagnostic combinations and assembles a barcode of desired length for a query taxon. First all single position DNCs are added, then based on the frequency of a nucleotide position in the DNCs of the 2 positions, and then based on the frequency of a position in longer DNCs
    OnePos = []
    positions_involved2 = []
    positions_involved_more = []
    for comb in Diagnostic_combinations:
        if len(comb) == 1:
            OnePos.append(comb[0])
        elif len(comb) == 2:
            for i in comb:
                positions_involved2.append(i)
        else:
            for j in comb:
                positions_involved_more.append(j)
    Setin = []
    for pos in sorted(positions_involved2, key=positions_involved2.count, reverse = True):
        if not pos in Setin:
            Setin.append(pos)
    for pos1 in sorted(positions_involved_more, key=positions_involved_more.count, reverse = True):
        if not pos1 in Setin:
            Setin.append(pos1)
    Ordered = OnePos + Setin
    return Ordered[:length]

def Screwed_dataset3(dataset, Taxon, sp_per_clade_to_screw):#implements a random_sequence function to build a random sequence dataset where 10% sequences of the dataset (but not more than 20 species per clade) are Pdiff percent different from the original one
    Pdiff_records=[]
    screwed_clades = {}
    seq_to_screw = random.sample(range(len(dataset)), len(dataset)/10)
    for i in range(len(dataset)):
        if i in seq_to_screw:
            clade_record = dataset[i][1]
            if clade_record not in list(screwed_clades.keys()):
                screwed_clades[clade_record] = 1                       
                Pdiff_record = [dataset[i][0], dataset[i][1], random_sequence(dataset[i][2], Percent_difference)]
                Pdiff_records.append(Pdiff_record)
            else:
                if screwed_clades[clade_record] < sp_per_clade_to_screw:
                    screwed_clades[clade_record] += 1                       
                    Pdiff_record = [dataset[i][0], dataset[i][1], random_sequence(dataset[i][2], Percent_difference)]
                    Pdiff_records.append(Pdiff_record)
                else:
                    Pdiff_records.append(dataset[i])
        else:
            Pdiff_records.append(dataset[i])
                
    Clades, clade_sorted_seqs, shared_positions, clade_sorted_spp = Step1(Pdiff_records)
    x,y,z = C_VP_PP(clade_sorted_seqs, Taxon, shared_positions, Cutoff)#STEP2
    return x, y

#3. IMPLEMENTATION###########################################################################################################################
#3.1. Analysis parameters
FragmentLen = 1041
"""
Fragment lengths:
COI - 658
CO2 - 684
ND2 - 1041
ND3 - 350
"""
Cutoff = 100
N1 = 10000
MaxLen1 = 12
Percent_difference = 2.5#function Screwed_dataset3
Sp_per_clade_to_screw = 10#function Screwed_dataset3


qCLADE = 'Synallaxis'#enter the clade to be examined
#3.2. READS MAIN DATA FILE and stores the sequences from the questioned clade in a separate list from remaining seqs
f = open("E:\Python\Molecular_diagnoses\Furnariidae_genera_ND2-GB_for_mold_corr.txt", "r")
imported=[]#set up a new dictionary with species and identifiers
for line in f:
    line=line.rstrip()#discards empty lines and Etc.
    words=line.split()
    if len(words) != 3 or len(words[2]) != FragmentLen or words[2].count('N')>5:
        continue
    else:
        imported.append(words) #indicates that first word in the line becomes an ID

Clades1, clade_sorted_seqs1, shared_positions1, clade_sorted_spp1 = Step1(imported)#Functions in 1.1. to compile full dataset for the Totat_Check
x1,y1,z1 = C_VP_PP(clade_sorted_seqs1, qCLADE, shared_positions1, Cutoff)#Functions in 1.2 to compile full dataset for the Totat_Check

All_clades = []
All_clades_species = {}
for i in imported:
    Clade_record = i[1]
    if Clade_record not in All_clades:
        All_clades.append(Clade_record)
        All_clades_species[Clade_record] = [i[0]]
    else:
       if i[0] in All_clades_species[Clade_record]:
           continue
       else:
           All_clades_species[Clade_record].append(i[0])
Totalspp = sum([len(i) for i in list(All_clades_species.values())])

class SortedDisplayDict(dict):#this is only to get a likable formatting of the barcode
   def __str__(self):
       return "[" + ", ".join("%r: %r" % (key, self[key]) for key in sorted(self)) + "]"

h = open("E:\Python\Molecular_diagnoses\Resampling_output\WALL_CLADES_resampling_Clade_%s_ND2-Barcode2.txt" %qCLADE, "w")
#h = open("E:\Python\Molecular_diagnoses\Resampling_output\WALL_CLADES_Synallaxis_qsp.txt", "w")
Sispecies = ['Synallaxis_propinqua', 'Schoeniophylax_phryganophilus', 'Certhiaxis_mustelinus', 'Certhiaxis_cinnamomeus']
print >>h, 'Cutoff', Cutoff
print >>h, 'N1', N1
print >>h, 'MaxLen1', MaxLen1
print >>h, 'Condition: len(Barcode_scores) >= 2 and Barcode_scores[-2] >= 90 and Barcode_scores[-1] >= 90'
print >>h, 'nsp', 'Nsp', 'len(ND_combinations)', 'N_involved_positions', 'sDNC', 'Valid?'
resampling_array = [i for i in range(2, 21, 2)]+[j for j in range(25, len(All_clades_species[qCLADE])+1, 5)]
for n in range(2, len(All_clades_species[qCLADE])+1, 2):#This block of code samples different number of species from all clades, this number is proportional to the number of species in a clade.
    N = 10
    while N > 0:
        qsp = []
        for Clade in All_clades:
            m = int(n*(float(len(All_clades_species[Clade]))) / len(All_clades_species[qCLADE])) + 1
            if len(All_clades_species[Clade]) > m:
                Csp = random.sample(All_clades_species[Clade], m)
                for i in Csp:
                    qsp.append(i)
            else:
                for j in All_clades_species[Clade]:
                    qsp.append(j)
        raw_records = []
        for k in imported:
            if k[0] in qsp:
                raw_records.append(k)
            else:
                continue
#        print n, len(qsp)
        for sisp in Sispecies:
            if not sisp in qsp:
                qsp.append(sisp)
        Clades, clade_sorted_seqs, shared_positions, clade_sorted_spp = Step1(raw_records)#Functions in 1.1.
        x,y,z = C_VP_PP(clade_sorted_seqs, qCLADE, shared_positions, Cutoff)#Functions in 1.2
        ND_combinations= []
        M = 5
        while M > 0:
            try:
                q = Diagnostic_combinations(qCLADE, x, y, N1, MaxLen1)#Function in 1.4
            except IndexError:
                print n, 'IndexError'
                continue
            for comb in q:
                if not comb in ND_combinations:
                    ND_combinations.append(comb)
            M-=1
        ND_combinations.sort(key=len)
        Allpos = []
        for comb in ND_combinations:
            for pos in comb:
                if not pos in Allpos:
                    Allpos.append(pos)
#        print >>h, n, N, len(ND_combinations), len(Allpos)
        print n, N, 'Number of retrieved diagnostic combinations:', len(ND_combinations), 'positions involved:', len(Allpos)

        Barcode_scores = []
        npos = len(ND_combinations[0])
        while npos <= min([25, len(Allpos)]):
            Barcode = GenerateBarcode2(ND_combinations, npos)
            Barcode_score = 0
            L = 100
            while L > 0:
                NComplist, NCPP = Screwed_dataset3(raw_records, qCLADE, Sp_per_clade_to_screw)
                NBarcode = [i for i in Barcode if i in list(NCPP.keys())]#i
                if ConditionD(NBarcode, NComplist, NCPP) == True:
                    Barcode_score +=1
                L -=1
            print npos, 'Barcode_score (100):', Barcode_score
            Barcode_scores.append(Barcode_score)
            if len(Barcode_scores) >= 3 and Barcode_scores[-2] >= 90 and Barcode_scores[-1] >= 90 and Barcode_scores[-3] >= 90:
                barcode1 = [i for i in Barcode if i in shared_positions1[qCLADE]]
                Valid = 0
                if ConditionD(barcode1, x1, y1) == True:
                    Valid += 1
                print SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in Barcode]}), Valid, '\n'
                print >>h, n, len(qsp), len(ND_combinations), len(Allpos), SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in Barcode]}), Valid
                break
            elif len(Barcode_scores) >= 10 and Barcode_scores[-5] > 70  and Barcode_scores[-4] >70 and Barcode_scores[-3] > 70 and Barcode_scores[-2] > 70 and Barcode_scores[-1] > 70:#'desperation barode' in the case of the key position
                barcode1 = [i for i in Barcode if i in shared_positions1[qCLADE]]
                Valid = 0
                if ConditionD(barcode1, x1, y1) == True:
                    Valid += 1
                print h, n, len(qsp), len(ND_combinations), len(Allpos), SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in Barcode]}), Valid
                print >>h, n, len(qsp), len(ND_combinations), len(Allpos), SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in Barcode]}), Valid
                break
            elif len(Barcode) == len(Allpos):
                barcode1 = [i for i in Barcode if i in shared_positions1[qCLADE]]
                Valid = 0
                if ConditionD(barcode1, x1, y1) == True:
                    Valid += 1
                print SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in Barcode]}), Valid, '\n'
                print >>h, n, len(qsp), len(ND_combinations), len(Allpos), SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in Barcode]}), Valid
                break
            else:
                npos += 1
        N-=1   
h.close()
f.close()