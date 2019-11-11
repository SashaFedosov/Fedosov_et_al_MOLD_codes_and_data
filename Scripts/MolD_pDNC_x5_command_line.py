"""
This script compiles pDNC-based DNA diagnoses for a pre-defined taxa in a dataset
This is working version to be used for both, the species-level and supraspecific taxa
"""
import os
import sys
import argparse
import random

#COMMON FUNCTIONS
#***STEP 1 - SORTING ENTRIES BY CLADE AND IDENTIFYING NUCLEOTIDE POSITIONS SHARED WITHIN CLADE
def Step1(raw_records):
    Clades=[]
    for i in range(len(raw_records)):
        Clade=raw_records[i][1]
        if Clade not in Clades:
            Clades.append(Clade)
    clade_sorted_seqs = {}
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
            if shared_nucleotide == True and csm != 'N':
                sh_pos.append(i)
        shared_positions[key]=sh_pos
    return Clades, clade_sorted_seqs, shared_positions

#***STEP 2 COMPILING COMPARISON LISTS FOR CLADES AND IDENTIFYING VARIABLE POSITIONS AND N PRIORITY POSITIONS WITH LARGEST CUTOFFS
def C_VP_PP(clade_sorted_seqs, clade, shared_positions, CUTOFF):# complist_variable_positions_priority_positions; Arguments: dictionary, string, dictionary
    CShN={}#a dictionary keys - clade shared positions, values - nucleotides at those positions
    for pos in shared_positions[clade]:
        CShN[pos] = clade_sorted_seqs[clade][0][pos]#creates a dictionary shared position : nucleotide
    complist=[]
    for key in clade_sorted_seqs:
        if key != clade:
            complist = complist + clade_sorted_seqs[key]#creates a list of all other sequences for comparison
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
    return complist, Clade_priority_positions, cutoffs

#***STEPS 3 RANDOM SEARCH ACROSS PRIORITY POSITIONS TO FIND RAW DIAGNOSTIC COMBINATIONS AND TO SUBSEQUENTLY REFINE THEM
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


def Diagnostic_combinations(qCLADE, complist, CPP, n1, maxlen1, maxlen2):#PUTS EVERYTHING TOGETHER - 20000 ROUNDS OF RANDOM SEARCH FOLLOWED BY REFINING OF 500 SHORTEST COMBINATIONS
    Achecked_ind = []
    bestlists = []
    n = n1
    while n>0:
        m = step_reduction_complist(qCLADE, complist, CPP, Achecked_ind)
        if len(m) < maxlen1 and sorted(m) not in bestlists:
            bestlists.append(sorted(m))
        n=n-1
    bestlists.sort(key=len)
    priority_lists = bestlists[:500]
    diagnostic_combinations = []
    for raw_comb in priority_lists:
        refined_comb = RemoveRedundantPositions(raw_comb, complist, CPP)
        if not refined_comb in diagnostic_combinations and len(refined_comb) <=maxlen2:
            diagnostic_combinations.append(refined_comb)
    diagnostic_combinations.sort(key=len)
    return diagnostic_combinations

###***STEP 4 ANALYSIS OF THE OUTPUT DNCs
def IndependentKey(diagnostic_combinations):
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
    independent_combinations.sort(key=len)
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

###########READ IN PARAMETER FILE AND DATA FILE
def get_args(): #arguments needed to give to this script
    parser = argparse.ArgumentParser(description="run MolD")
    required = parser.add_argument_group("required arguments")
    required.add_argument("-i", help="textfile with parameters of the analysis", required=True)
    return parser.parse_args()

def main():
    args = get_args()
    ParDict = {}
    with open(args.i) as params:
        for line in params:
            line = line.strip()
            if line.startswith('#'):
                pass
            else:
                if len(line.split('=')) == 2 and len(line.split('=')[1]) != 0:
                    ParDict[line.split('=')[0]] = line.split('=')[1]

    f = open(ParDict['INPUT_FILE'], 'r') 
    imported=[]#set up a new dictionary with species and identifiers
    for line in f:
        line=line.rstrip()
        words=line.split()
        if len(words) != 3:
            print 'Check number of entries in', words[0]
            #break
        else:
            imported.append([words[0], words[1], words[2].upper()])
    f.close()
    if len(set([len(i[2]) for i in imported])) != 1:
        print 'Alignment contains sequences of different lengths:', set([len(i[2]) for i in imported])
    else:
        FragmentLen = len(imported[0][2])
    if 'NumberN' in list(ParDict.keys()):#How many ambiguously called nucleotides are allowed
        NumberN = int(ParDict['NumberN'])
    else:
        NumberN = 5
    raw_records=[]
    for i in imported:
        if i[2].count('N') < NumberN and len(i[2]) == FragmentLen:
            raw_records.append(i)
    print 'Maximum undetermined nucleotides allowed:', NumberN
    print 'Length of the alignment:', FragmentLen    
    print 'Read in', len(raw_records), 'sequences'
   
    #############################################READ IN OTHER ANALYSIS PARAMETERS
    if ParDict['qTAXA'] in ['ALL', 'All', 'all']:#qTAXA
        qCLADEs = []
        for i in raw_records:
            if not i[1] in qCLADEs:
                qCLADEs.append(i[1])
    elif ParDict['qTAXA'][0] == '>':
        NumSeq = int(ParDict['qTAXA'][1:])
        Taxarecords = [i[1] for i in raw_records]
        qCLADEs = []
        for j in Taxarecords:
            if Taxarecords.count(j) >= NumSeq and not j in qCLADEs:
                qCLADEs.append(j)
    else:
        qCLADEs = ParDict['qTAXA'].split(',')
    print 'focus taxa:', qCLADEs, len(qCLADEs)

    if 'Cutoff' in list(ParDict.keys()):#CUTOFF Number of the informative positions to be considered, default 100
        Cutoff = int(ParDict['Cutoff'])
    else:
        Cutoff = 100
    print 'Cutoff set as:', Cutoff
    if 'Number_of_iterations' in list(ParDict.keys()):#Number iterations of MolD
        N1 = int(ParDict['Number_of_iterations'])
    else:
        N1 = 10000
    print 'Number iterations of MolD set as:', N1
    
    if 'MaxLen1' in list(ParDict.keys()):#Maximum length for the raw pDNCs
        MaxLen1 = int(ParDict['MaxLen1'])
    else:
        MaxLen1 = 12
    print 'Maximum length of raw pDNCs set as:', MaxLen1
    
    if 'MaxLen2' in list(ParDict.keys()):#Maximum length for the refined pDNCs
        MaxLen2 = int(ParDict['MaxLen2'])
    else:
        MaxLen2 = 7
    print 'Maximum length of refined pDNCs set as:', MaxLen2

    ###################################################IMPLEMENTATION
    class SortedDisplayDict(dict):
        def __str__(self):
            return "[" + ", ".join("%r: %r" % (key, self[key]) for key in sorted(self)) + "]"
    Clades, clade_sorted_seqs, shared_positions = Step1(raw_records)#STEP1
    for qCLADE in qCLADEs:
        print '\n************', qCLADE, '************'
        print 'Sequences analyzed:', len(clade_sorted_seqs[qCLADE])
        x,y,z = C_VP_PP(clade_sorted_seqs, qCLADE, shared_positions, Cutoff)#STEP2
        ND_combinations= []
        N = 5
        while N > 0: #STEP3 repeated 5 times   
            try:
                q = Diagnostic_combinations(qCLADE, x, y, N1, MaxLen1, MaxLen2)
            except IndexError:
                print N, 'IndexError'
                continue
            for comb in q:
                if not comb in ND_combinations:
                    ND_combinations.append(comb)
            N-=1
            print N
        ND_combinations.sort(key=len)
        try:
            r, s = IndependentKey(ND_combinations)
        except IndexError:
            print 'Sorry, not a single pDNC found for', qCLADE
            continue
        print 'Independent pDNCs:', len(r)
        for comb in r:
            print SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in comb]})
        print 'Key positions:', len(s), SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in s]})
        print '25 shortest pDNCs:'
        for comb in ND_combinations[:min([25, len(ND_combinations)-1])]:
            print SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in comb]})

if __name__ == "__main__":
	main()