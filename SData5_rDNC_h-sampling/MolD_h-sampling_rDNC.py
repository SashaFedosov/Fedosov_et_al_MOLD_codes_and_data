"""
This is the script used for the sDNC resampling of the query haplotypes (species composition of the partial datasets remains unchanged) (second group of resampling results, figs 6A - C)
"""
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

#1.2. COMPILING COMPARISON LISTS FOR CLADES AND IDENTIFYING VARIABLE POSITIONS AND N PRIORITY POSITIONS WITH LARGEST CUTOFFS
def C_VP_PP(clade_sorted_seqs, clade, shared_positions, CUTOFF):# complist_variable_positions_priority_positions; Arguments: dictionary, string, dictionary
    CShN={}#a dictionary keys - clade shared positions, values - nucleotides at those positions
    for pos in shared_positions[clade]:
        CShN[pos] = clade_sorted_seqs[clade][0][pos]#creates a dictionary shared position : nucleotide
    complist=[]
    for key in clade_sorted_seqs:
        if key != clade:
            complist = complist + clade_sorted_seqs[key]#creates a list of all other sequences for comparison
    cutoffs = {}
    pures = []####! newline
    for key in CShN:
        newcomplist = []
        for k in complist:
            if k[key] == CShN[key]:
                newcomplist.append(k)
            else: continue
        cutoffs[key] = len(complist) - len(newcomplist)
        if len(newcomplist) == 0:####! newline
            pures.append(key)####! newline
    CPP = []
    for key in sorted(cutoffs, key = cutoffs.get, reverse = True):
        CPP.append(key)
    if str(CUTOFF)[0] == '>':#VERYNEW
        Clade_priority_positions = {pos:CShN[pos] for pos in CPP if cutoffs[pos] > int(CUTOFF[1:])}#VERYNEW
    else:#VERYNEW
        Clade_priority_positions = {}
        for position in CPP[:int(CUTOFF)]:#Here you define how many of the clade shared combinations are used in subsequent search
            Clade_priority_positions[position] = CShN[position]
    return complist, Clade_priority_positions, cutoffs, pures####! pures added

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
def Diagnostic_combinations(qCLADE, complist, CPP, n1, maxlen1, maxlen2):
    Achecked_ind = []
    bestlists = []
    n = n1
    while n>0:#STEP3 proposes raw diagnostic combinations
        m = step_reduction_complist(qCLADE, complist, CPP, Achecked_ind)
        if len(m) < maxlen1 and sorted(m) not in bestlists:
            bestlists.append(sorted(m))
        n=n-1
    bestlists.sort(key=len)
    priority_lists = bestlists#[:500]
    diagnostic_combinations = []
    for raw_comb in priority_lists:
        refined_comb = RemoveRedundantPositions(raw_comb, complist, CPP)
        if not refined_comb in diagnostic_combinations and len(refined_comb) <=maxlen2:
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

#1.5 FUNCTIONS TO BUILD AND REFINE sDNC
def PositionArrays(Motifs):#VERYNEW
    PositionArrays = []
    VarPosList = []
    for i in range(len(Motifs[0])):
        Const = True
        array = [Motifs[0][i]]
        for j in range(len(Motifs[1:])):
            if Motifs[j][i] != 'N':
                if Motifs[j][i] != array[-1]:
                    Const = False
                array.append(Motifs[j][i])
        PositionArrays.append(array)
        if Const == False:
            VarPosList.append(i)
    return PositionArrays, VarPosList

def random_sequence_new(SEQ, PositionArrays, VarPosList, Pdiff):#VERYNEW
    n = len(SEQ)*Pdiff/100
    N = random.sample(range(1, n), 1)[0]
    PosToChange = random.sample(VarPosList, N)
    NEWSEQ = ''
    for i in range(len(SEQ)):
        if i not in PosToChange:
            NEWSEQ = NEWSEQ + SEQ[i]
        else:
            newarray = [j for j in PositionArrays[i] if j != SEQ[i]]
            newbase = random.sample(newarray, 1)[0]
            NEWSEQ = NEWSEQ + newbase
    return NEWSEQ

def GenerateBarcode_new(Diagnostic_combinations, length):#VERYNEW #This function calculates diagnostic combinations and assembles a barcode of desired length for a query taxon. First all single position DNCs are added, then based on the frequency of a nucleotide position in the DNCs of the 2 positions, and then based on the frequency of a position in longer DNCs
    len1 = []
    len2 = []
    lenmore = []
    for comb in Diagnostic_combinations:
        if len(comb) == len(Diagnostic_combinations[0]):
            for i in comb:
                len1.append(i)
        elif len(comb) == len(Diagnostic_combinations[0])+1:
            for j in comb:
                len2.append(j)
        else:
            for k in comb:
                lenmore.append(k)
    if len(Diagnostic_combinations[0]) == 1:
        Setin = len1
    else:
        Setin = []
        for pos in sorted(len1, key=len1.count, reverse = True):
            if not pos in Setin:
                Setin.append(pos)
    for pos1 in sorted(len2, key=len2.count, reverse = True):
        if not pos1 in Setin:
            Setin.append(pos1)
    for pos2 in sorted(lenmore, key=lenmore.count, reverse = True):
        if not pos2 in Setin:
            Setin.append(pos2)
    return Setin[:length]

def GenerateBarcode_new2(Diagnostic_combinations, length):#VERYNEW #This function calculates diagnostic combinations and assembles a barcode of desired length for a query taxon. First all single position DNCs are added, then based on the frequency of a nucleotide position in the DNCs of the 2 positions, and then based on the frequency of a position in longer DNCs
    Setin = Diagnostic_combinations[0]
    if len(Setin) < length:
        for comb in Diagnostic_combinations[1:]:
            for pos in comb:
                if pos in Setin:
                    continue
                else:
                    if len(Setin) < length:
                        Setin.append(pos)
                    else:
                        return Setin
    else:
        return Setin

def Screwed_dataset_new(raw_records, nseq_per_clade_to_screw, PositionArrays, VarPosList, Percent_difference, Taxon, Cutoff):#VERYNEW
    clades=[]
    for i in range(len(raw_records)):
        Clade=raw_records[i][1]
        if Clade not in clades:
            clades.append(Clade)
    clade_sorted_seqs = {}
    for letter in clades:
        clade_sorted_seqs[letter]=[]
        for i in range(len(raw_records)):
            if raw_records[i][1]==letter:
                clade_sorted_seqs[letter].append(raw_records[i][2])
    
    for clade in clades:
        seqlist = clade_sorted_seqs[clade]
        newseqs = []
        if len(seqlist) > nseq_per_clade_to_screw:
            iSTS = random.sample(range(len(seqlist)), nseq_per_clade_to_screw)
            for k in range(len(seqlist)):
                if k in iSTS:
                    newseq = random_sequence_new(seqlist[k], PositionArrays, VarPosList, Percent_difference)
                else:
                    newseq = seqlist[k]
                newseqs.append(newseq)
        elif len(clade_sorted_seqs[clade]) == nseq_per_clade_to_screw:
            for k in range(len(seqlist)):
                newseq = random_sequence_new(seqlist[k], PositionArrays, VarPosList, Percent_difference)
                newseqs.append(newseq)
        else:
            for i in range(nseq_per_clade_to_screw):
                seq = random.sample(seqlist, 1)[0]
                newseq = random_sequence_new(seq, PositionArrays, VarPosList, Percent_difference)
                newseqs.append(newseq)
        clade_sorted_seqs[clade] = newseqs
    
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
    
    x,y,z,pures = C_VP_PP(clade_sorted_seqs, Taxon, shared_positions, Cutoff)#STEP2####
    return x, y

def isSameHap(seq1, seq2):
    issame = True
    for k in range(len(seq1)):
        if seq1[k] == seq2[k] or seq1[k] == 'N' or seq2[k] == 'N':
            continue
        else:
            issame = False
            break
    return issame

def Hcount(records):#this function returns number of the unique haplotype in the alignment, ignoring 'N's
    Haps = [records[0]]
    for thing in records[1:]:
        new = True
        for hap in Haps:
            if isSameHap(thing, hap) == True:
                if thing.count('N') < hap.count('N'):
                    hap = thing
                new = False
                break
        if new == True:
            Haps.append(thing)
    return Haps

def DefineVarPos(seqs):
    VarPos = []
    for i in range(len(seqs[0])):
        nlist = [seq[i] for seq in seqs if seq[i] != 'N']
        if len(nlist) != 0:
            if len(nlist) == nlist.count(nlist[0]): 
                continue
            else:
                VarPos.append(i)
    return VarPos

#3. IMPLEMENTATION###########################################################################################################################
#3.1. Analysis parameters

Cutoff = '>1'
N1 = 20000
MaxLen1 = 12
MaxLen2 = 5
Seq_per_clade_to_screw = 10#function Screwed_dataset_new
Percent_difference = 2#function Screwed_dataset_new
NumberN=50
threshold = 75

"""
Indespensables = {}
Indespensables['similis'] = ['sinensis']
Indespensables['laevis'] = ['cucullata']
Indespensables['melanica'] = ['pulex0','middendorffiana','cf.pulicaria', 'tenebrosa']
Indespensables['pulex0'] = ['melanica','middendorffiana','cf.pulicaria', 'tenebrosa']
Indespensables['miliaris'] = ['chaldaeus', 'ebraeus', 'aristophanes']
Indespensables['chaldaeus'] = ['judaeus', 'coronatus', 'ebraeus']
Indespensables['ebraeus'] = ['judaeus', 'coronatus', 'chaldaeus']
"""
#3.2. READS MAIN DATA FILE and stores the sequences from the questioned clade in a separate list from remaining seqs
f = open("D:\Daphnia_COI_657corr.txt", "r")
#f = open("D:\Conus_4sppCOI.txt", "r")
#f = open("D:\Xenuroturris_COI.txt", "r")

#qCLADE = 'longispina'#enter the clade to be examined
qCLADE = 'pulex'
gene = 'COI'

imported=[]#set up a new dictionary with species and identifiers
for line in f:
    line=line.rstrip()#discards empty lines and Etc.
    words=line.split()
    if len(words) != 3 or words[2].count('N') + words[2].count('n') + words[2].count('-') > NumberN:
        continue
    else:
        imported.append([words[0], words[1], words[2].upper().replace('-', 'D').replace('R', 'N').replace('Y', 'N').replace('S','N').replace('W','N').replace('M','N').replace('K','N')])
FragmentLen = len(imported[0][2])
Clades1, clade_sorted_seqs1, shared_positions1 = Step1(imported)#Functions in 1.1. to compile full dataset for the Totat_Check

UHaplotypes = {}
for clade in Clades1:
    UHaplotypes[clade] = Hcount(clade_sorted_seqs1[clade])
    #print clade, len(clade_sorted_seqs1[clade]), len(UHaplotypes[clade])

x1,y1,z1,pures1 = C_VP_PP(clade_sorted_seqs1, qCLADE, shared_positions1, Cutoff)#Functions in 1.2 to compile full dataset for the Totat_Check

class SortedDisplayDict(dict):#this is only to get a likable formatting of the barcode
   def __str__(self):
       return "[" + ", ".join("%r: %r" % (key, self[key]) for key in sorted(self)) + "]"
"""
h = open("D:\sDNC_qHap_resampling\\"+qCLADE+'_'+gene+"pdiff05.txt", "w")
print >>h, 'Cutoff', Cutoff
print >>h, 'N1', N1
print >>h, 'MaxLen1', MaxLen1
print >>h, 'Condition: len(Barcode_scores) >= 2 and Barcode_scores[-2] >= 75 and Barcode_scores[-1] >= 75\n'
print 'n', 'N', 'NSp', 'len(raw_records)', 'qPVS', 'dPVS', 'nPureDNCs', 'Retrieved_pDNCs', 'N_positions_in_all pDNCs'
print >>h, 'n', 'N', 'NSp', 'len(raw_records)', 'qPVS', 'dPVS', 'nPureDNCs', 'Retrieved_pDNCs', 'N_positions_in_all pDNCs', 'sDNC', 'Valid?'
"""
#for n in range(2, len(UHaplotypes[qCLADE])+1):#This block of code samples different number of species from all clades, this number is proportional to the number of species in a clade.
n = len(UHaplotypes[qCLADE])+1#This block of code samples different number of species from all clades, this number is proportional to the number of species in a clade.
if n != 0:
    N = 1
    while N > 0:
        dataset = open('D:\sDNC_qHap_resampling\\dataset_'+qCLADE+'_'+gene+'_'+str(n)+str(N)+'.txt', 'w')
        raw_records = []
        mnull = []
        for Clade in Clades1:
            m = int(n*(float(len(UHaplotypes[Clade])+0.4)) / len(UHaplotypes[qCLADE]))
            if m == 0:
                mnull.append(Clade)
            else:
                if len(UHaplotypes[Clade]) > m:
                    CHaps = random.sample(UHaplotypes[Clade], m)
                else:
                    CHaps = UHaplotypes[Clade]
                for thing in CHaps:
                    raw_records.append([Clade, Clade, thing])
                    print >>dataset, Clade, Clade, thing
        print >>dataset, '##########################################'
        for clade in mnull:
            Hap = random.sample(UHaplotypes[clade], 1)[0]
            raw_records.append([clade, clade, Hap])
            print >>dataset, clade, clade, Hap
        dataset.close()
        N -=1
"""        
        PosArrays, VarPosList = PositionArrays([i[2] for i in raw_records])
        Clades, clade_sorted_seqs, shared_positions = Step1(raw_records)#Functions in 1.1.
        SubSampleHaps = {}
        for smallclade in Clades:
            SubSampleHaps[smallclade] = Hcount(clade_sorted_seqs[smallclade])
        print n, N, len(set([i[1] for i in raw_records])), len(raw_records), float(len(DefineVarPos(SubSampleHaps[qCLADE])))/len(raw_records[0][2]), float(len(DefineVarPos([i[2] for i in raw_records])))/len(raw_records[0][2]),
        print >>h, n, N, len(set([i[1] for i in raw_records])), len(raw_records), float(len(DefineVarPos(SubSampleHaps[qCLADE])))/len(raw_records[0][2]), float(len(DefineVarPos([i[2] for i in raw_records])))/len(raw_records[0][2]),
        
        x,y,z,pures = C_VP_PP(clade_sorted_seqs, qCLADE, shared_positions, Cutoff)#Functions in 1.2
        newy = {key:y[key] for key in y if not key in pures}
        ND_combinations = [[item] for item in pures]
        print len(pures),
        if len(pures) < 10:
            try:
                q = Diagnostic_combinations(qCLADE, x, newy, N1, MaxLen1, MaxLen2)#Function in 1.4
            except IndexError:
                print n, 'IndexError'
                continue
            for comb in q:
                if not comb in ND_combinations:
                    ND_combinations.append(comb)
        ND_combinations.sort(key=len)
        
        Allpos = []
        for comb in ND_combinations:
            for pos in comb:
                if not pos in Allpos:
                    Allpos.append(pos)
        print len(ND_combinations), len(Allpos)
        print >>h, len(ND_combinations), len(Allpos),
        
        Barcode_scores = []
        npos = len(ND_combinations[0])
        BestBarcode = 'none'####! newline
        while npos <= min([10, len(Allpos)]):
            Barcode = GenerateBarcode_new(ND_combinations, npos)
            Barcode_score = 0
            L = 100
            while L > 0:
                NComplist, NCPP = Screwed_dataset_new(raw_records, Seq_per_clade_to_screw, PosArrays, VarPosList, Percent_difference, qCLADE, Cutoff)
                NBarcode = [i for i in Barcode if i in list(NCPP.keys())]#i
                if len(Barcode) - len(NBarcode) <= 1 and ConditionD(NBarcode, NComplist, NCPP) == True:####! new condition
                    Barcode_score +=1
                L -=1
            print npos, 'Barcode_score (100):', Barcode_score
            if Barcode_score >= threshold and len(Barcode_scores) >= 1 and Barcode_score >= Barcode_scores[-1]:
                BestBarcode = Barcode####!newline
            Barcode_scores.append(Barcode_score)    
            if len(Barcode_scores) >= 3 and Barcode_scores[-1] >= threshold and Barcode_scores[-2] >= threshold:
                barcode1 = [i for i in Barcode if i in shared_positions1[qCLADE]]
                Valid = 0
                if len(Barcode) - len(barcode1) <=1 and ConditionD(barcode1, x1, y1) == True:
                    Valid += 1
                print SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in Barcode]}), Valid, '\n'
                print >>h, n, len(raw_records), len(ND_combinations), len(Allpos), SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in Barcode]}), Valid
                break
            elif len(Barcode_scores) == min([10, len(Allpos)]) and BestBarcode != 'none':#'desperation barode' in the case of the key position
                barcode1 = [i for i in BestBarcode if i in shared_positions1[qCLADE]]
                Valid = 0
                if len(BestBarcode) - len(barcode1) <=1 and ConditionD(barcode1, x1, y1) == True:
                    Valid += 1
                print h, n, len(raw_records), len(ND_combinations), len(Allpos), SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in BestBarcode]}), Valid
                print >>h, n, len(raw_records), len(ND_combinations), len(Allpos), SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in BestBarcode]}), Valid
                break
            elif len(Barcode_scores) == min([10, len(Allpos)]) and BestBarcode == 'none':
                barcode1 = [i for i in BestBarcode if i in shared_positions1[qCLADE]]
                Valid = 0
                print SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in Barcode]}), Valid, '\n'
                print >>h, n, len(raw_records), len(ND_combinations), len(Allpos), SortedDisplayDict({pos : y[pos-1] for pos in [i+1 for i in Barcode]}), Valid
                break
            else:
                npos += 1
        N-=1   
h.close()
"""
f.close()
