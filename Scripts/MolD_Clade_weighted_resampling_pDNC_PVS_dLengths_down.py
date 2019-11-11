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
    clade_sorted_seqs={}
    for letter in Clades:
        clade_sorted_seqs[letter]=[]
        for i in range(len(raw_records)):
            if raw_records[i][1]==letter:
                clade_sorted_seqs[letter].append(raw_records[i][2])
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
    return Clades, clade_sorted_seqs, shared_positions

#1.2. COMPILING COMPARISON LISTS FOR CLADES AND IDENTIFYING VARIABLE POSITIONS AND N PRIORITY POSITIONS WITH LARGEST CUTOFFS
def C_VP_PP(clade_sorted_seqs, clade, shared_positions, CUTOFF):# complist_variable_positions_priority_positions; Arguments: dictionary, string, dictionary
    CShN={}
    for pos in shared_positions[clade]:
        CShN[pos] = clade_sorted_seqs[clade][0][pos]#creates a dictionary shared position : nucleotide
    complist=[]
    for key in clade_sorted_seqs:
        if key != clade:
            complist = complist + clade_sorted_seqs[key]#creates a list of all other sequences for comparison
    cstrings={}
    for key in CShN:
        cstrings[key] = ''
        for j in complist:
            cstrings[key] = cstrings[key] + j[key]
    variable_positions = {}
    for key in cstrings:
        if cstrings[key].count(CShN[key]) < len(cstrings[key]):
            variable_positions[key] = cstrings[key]
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

#1.5 AFUNCTION TO VERIFY WHETHER THE RETURNED DIAGNOSIS COMPUTED FOR THE SUBSET OF TAXA HOLDS WHEN TESTED FOR COMPLETE SET OF TAXA IN A CLADE
def Total_Check(diagnosis, complist, CShN):#list of lists, list, dictionary
    WORKS1 = 0
    WORKS2 = 0
    WORKS3 = 0
    WORKS4 = 0
    WORKSMORE = 0
    for comb in diagnosis:
        Good_comb = True
        for pos in comb:
            if pos in list(CShN.keys()):
                continue
            else:
                Good_comb = False
                break
        if Good_comb == True and ConditionD(comb, complist, CShN) == True and len(comb) == 1:
            WORKS1+=1
        elif Good_comb == True and ConditionD(comb, complist, CShN) == True and len(comb) == 2:
            WORKS2+=1
        elif Good_comb == True and ConditionD(comb, complist, CShN) == True and len(comb) == 3:
            WORKS3+=1
        elif Good_comb == True and ConditionD(comb, complist, CShN) == True and len(comb) == 4:
            WORKS4+=1
        elif Good_comb == True and ConditionD(comb, complist, CShN) == True and len(comb) > 4:
            WORKSMORE+=1
        else:
            continue
    return [WORKS1, WORKS2, WORKS3, WORKS4, WORKSMORE]

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
qCLADE = 'Synallaxis'#enter the clade to be examined
#3.2. READS MAIN DATA FILE and stores the sequences from the questioned clade in a separate list from remaining seqs
f = open("E:\Python\Molecular_diagnoses\Furnariidae_genera_ND2_for_mold_corr.txt", "r")
imported=[]#set up a new dictionary with species and identifiers
for line in f:
    line=line.rstrip()#discards empty lines and Etc.
    words=line.split()
    if len(words) != 3 or len(words[2]) != FragmentLen or words[2].count('N')>5:
        continue
    else:
        imported.append(words) #indicates that first word in the line becomes an ID
"""
#This block of code is executed when the full dataset comes from the separate file
g = open("E:\Python\Molecular_diagnoses\Furnariidae_SF_ND2-GB_for_mold_corr.txt", "r")
imported1=[]#set up a new dictionary with species and identifiers
for line in g:
    line=line.rstrip()#discards empty lines and Etc.
    words=line.split()
    if len(words) != 3 or len(words[2]) != FragmentLen or words[2].count('N')>5:
        continue
    else:
        imported1.append(words) #indicates that first word in the line becomes an ID
"""
Clades1, clade_sorted_seqs1, shared_positions1 = Step1(imported)#Functions in 1.1. to compile full dataset for the Totat_Check
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


h = open("E:\Python\Molecular_diagnoses\Resampling_output\WALL_CLADES_resampling_PVS_Clade_%s_downsampled_COI_cutoff_100_5x10k_lengths.txt" %qCLADE, "w")
print >>h, 'Cutoff', Cutoff
print >>h, 'N1', N1
print >>h, 'MaxLen1', MaxLen1
print >>h, 'Nsp', 'D_var_sites', 'Len_distribution', 'WORKS_len_distribution', 'Shortest_combination', 'n_Independent_combinations', 'n_Key_positions'
#resampling_array = [i for i in range(2, 51)]+[j for j in range(55, len(All_clades_species[qCLADE])+1, 5)]
for n in range(2, len(All_clades_species[qCLADE])+1):#This block of code samples different number of species from all clades, this number is proportional to the number of species in a clade.
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
        sp_in = []#for a reduced sampling
        for k in imported:
            if k[0] in qsp and k[0] not in sp_in:#for a reduced sampling
                raw_records.append(k)
                sp_in.append(k[0])#for a reduced sampling
            else:
                continue
        Clades, clade_sorted_seqs, shared_positions = Step1(raw_records)#Functions in 1.1.
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
        try:
            r, s = IndependentKey(ND_combinations)
        except IndexError:
            print 'Sorry, not a single diagnostic nucleotide combination found'
            continue
        diagnosis = [0,0,0,0,0]
        for comb in ND_combinations:
            if len(comb) < 5:
                diagnosis[len(comb)-1]+=1
            else:
                diagnosis[4]+=1

        print >>h, n, str(float(FragmentLen - len(shared_positions[qCLADE])) / FragmentLen), diagnosis, Total_Check(ND_combinations, x1, z1), len(ND_combinations[0]), len(r), len(s)
        print n, N, len(qsp), str(float(FragmentLen - len(shared_positions[qCLADE])) / FragmentLen), diagnosis, '(', Total_Check(ND_combinations, x1, z1), ')'
        N-=1
h.close()
f.close()