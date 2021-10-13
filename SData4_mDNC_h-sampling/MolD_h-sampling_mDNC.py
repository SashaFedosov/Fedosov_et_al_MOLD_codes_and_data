import os, sys
import random
"""
This script performs haplotype resampling for a selected query taxon, and checkes the validity of each of the resulting pDNC against the final dataset.
The pDNCs are checked in several bins, corresponding to lengths of DNC (1-5 positions), and whether or not DNC contains AA-defining position(s)
Taxa resampling: proportional to the number of haplotypes per species, but all taxa included right away.
Status:  This is the final version. This script produced all pDNC resampling results (figs 3,4).
Output:
    1.Each generated individual dataset is recorded into a file (available upon request)
    2.Main output named by species name in the q taxa pDNC_h-resampling directory is plotted by one of the two plotting scripts to produce
    curves of changing pDNC reliability and reliability scatterplots.
"""
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
        for position in CPP[:CUTOFF]:#Here you define how many of the clade shared combinations are used in subsequent search
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

def ImportTrTab(filename):#VERYNEW
    data = open(filename, 'r')
    AAs = data.readline().strip().replace(' ', '').split('=')[1]
    first = data.readline().strip().replace(' ', '').split('=')[1]
    second = data.readline().strip().replace(' ', '').split('=')[1]
    third = data.readline().strip().replace(' ', '').split('=')[1]
    data.close()
    TrTable = {}
    for i in range(len(first)):
        Codon = first[i]+second[i]+third[i]
        TrTable[Codon] = AAs[i]
    return TrTable

def TranslateNSeq(Seq, start, TrTable):#VERYNEW
    AAseq = ''
    for i in range(start, len(Seq)-2, 3):
        codon = Seq[i:i+3]
        if codon[0] in 'ACGT' and codon[1] in 'ACGT' and codon[2] in 'ACGT':
            AAseq = AAseq + TrTable[codon]
        else:
            AAseq = AAseq + 'X'
    return AAseq 

def RevCompl(seq):#VERYNEW
    reseq = ''
    complems = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'R':'N', 'Y':'N', 'S':'N', 'W':'N', 'M':'N', 'K':'N'}
    for i in seq[::-1]:
        reseq = reseq+complems[i]
    return reseq

def IdentifyStart(Seqs, TrTable):
    BestNStops = len(Seqs[0])*len(Seqs)
    CodonP = 1
    RC = False
    for start in [0, 1, 2]:
        Nstops = 0
        for seq in Seqs:
            Nstops += TranslateNSeq(seq, start, TrTable).count('*')
        if Nstops < BestNStops:
            BestNStops = Nstops
            CodonP = start + 1
    for start in [0, 1, 2]:
        Nstops = 0
        for seq in Seqs:
            Nstops += TranslateNSeq(RevCompl(seq), start, TrTable).count('*')
        if Nstops < BestNStops:
            BestNStops = Nstops
            CodonP = 0 - start - 1
            RC = True
    return CodonP, RC, BestNStops
            
def CodonPositionLists(lenseq, First_index):#VERYNEW
    Fst = range(lenseq)[First_index::3]
    if Fst[0] == 2:
        Snd = range(lenseq)[::3]
    else:    
        Snd = range(lenseq)[First_index+1::3]
    if Snd[0] == 2:
        Trd = range(lenseq)[::3]
    else:
        Trd = range(lenseq)[First_index+2::3]
    return Fst, Snd, Trd

def ExtractUniqueCodons(pos, Seqs, Fst, Snd, Trd):
    seqlen = len(Seqs[0])
    codons = []
    if pos in Fst and pos in range(seqlen-2):
        for seq in Seqs:
            codon = seq[pos:pos+3]
            codons.append(codon)
    elif pos in Snd and pos in range(1, seqlen-1):
        for seq in Seqs:
            codon = seq[pos-1:pos+2]
            codons.append(codon)
    elif pos in Trd and pos in range(2, seqlen):
        for seq in Seqs:
            codon = seq[pos-2:pos+1]
            codons.append(codon)
    else:
        codon = 'NNN'
        codons.append(codon)
    Ucodons = []
    for cod in codons:
        if not cod in Ucodons:
            Ucodons.append(cod)
    return Ucodons

def SynNonSyn(codons, TrTable):#VERYNEW
    noNcodons = [codon for codon in codons if not 'N' in codon]
    syn = True
    if len(noNcodons) >= 2:
        AA = TrTable[noNcodons[0]]
        for i in range(1, len(noNcodons)):
            if TrTable[noNcodons[i]] != AA and TrTable[noNcodons[i]] != '*':
                syn = False
                break
    return syn

def RatePositionNucl(pos, Seqs, Fst, Snd, Trd, TrTable):#VERYNEW REPLACE
    codons = ExtractUniqueCodons(pos, Seqs, Fst, Snd, Trd)
    InvTrTable = {}
    for thing in TrTable.values():
        if not thing in InvTrTable.keys():
            InvTrTable[thing] = [elem for elem in TrTable.keys() if TrTable[elem] == thing]
    rate = []
    for codon in codons:
        if codon not in TrTable.keys():
            continue
        else:
            AAseq = TrTable[codon]
            if len(InvTrTable[AAseq]) > 1:
                for thing in InvTrTable[AAseq]:
                    if pos in Fst and codon[1] == thing[1] and codon[2] == thing[2] and thing != codon and not thing[0] in rate:
                        rate.append(thing[0])
                    elif pos in Snd and codon[0] == thing[0] and codon[2] == thing[2] and thing != codon and not thing[1] in rate:
                        rate.append(thing[1])
                    elif pos in Trd and codon[0] == thing[0] and codon[1] == thing[1] and thing != codon and not thing[2] in rate:
                        rate.append(thing[2])
                    else:
                        continue
    return rate

def PosRatesDictNucl(CPP, Seqs, Fst, Snd, Trd, TrTable):#VERYNEW
    RateDict = {}
    for key in CPP:
        if key != 0 and key != 1:
            RateDict[key] = RatePositionNucl(key, Seqs, Fst, Snd, Trd, TrTable)
        else:
            RateDict[key] = [i for i in 'ACGT' if i != CPP[key]]
    return RateDict

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
NumberN=50
qCLADE = 'longispina'#enter the clade to be examined
#qCLADE = sys.argv[1]
#3.2. READS MAIN DATA FILE and stores the sequences from the questioned clade in a separate list from remaining seqs
#f = open("D:\Xenuroturris_COI.txt", "r")
f = open("D:\Daphnia_COI_657corr.txt", "r")
#f = open("D:\Conus_4sppCOI.txt", "r")
#f = open("D:\Tanytarsus_PGDI.txt", "r")
gene = 'COI'
imported=[]#set up a new dictionary with species and identifiers
for line in f:
    line=line.rstrip()#discards empty lines and Etc.
    words=line.split()
    if len(words) != 3 or words[2].count('N') + words[2].count('n') + words[2].count('-') > NumberN:
        continue
    else:
        imported.append([words[0], words[1], words[2].upper().replace('-', 'N')])#NEW

PosArrays, VarPosList = PositionArrays([i[2] for i in imported])###VERYNEW
FragmentLen = len(imported[0][2])

CurrTrTable = 'D:\TrTable5.txt'
TrTable = ImportTrTab(CurrTrTable)

# lET'S GO
First_index = IdentifyStart([i[2] for i in imported], TrTable)[0] -1
FCP, SCP, TCP = CodonPositionLists(FragmentLen, First_index)

Clades1, clade_sorted_seqs1, shared_positions1 = Step1(imported)#Functions in 1.1. to compile full dataset for the Totat_Check
AllHaps = 0
UHaplotypes = {}
for clade in Clades1:
    UHaplotypes[clade] = Hcount(clade_sorted_seqs1[clade])
    AllHaps += len(UHaplotypes[clade])
    #print clade, len(clade_sorted_seqs1[clade]), len(UHaplotypes[clade])

x1,y1,z1,pures1 = C_VP_PP(clade_sorted_seqs1, qCLADE, shared_positions1, len(shared_positions1[qCLADE]))#Functions in 1.2 to compile full dataset for the Totat_Check

class SortedDisplayDict(dict):#this is only to get a likable formatting of the barcode
   def __str__(self):
       return "[" + ", ".join("%r: %r" % (key, self[key]) for key in sorted(self)) + "]"

h = open("D:\qHaps_pDNC_resampling\\"+qCLADE+'_'+gene+".txt", "w")
print >>h, 'Cutoff', Cutoff
print >>h, 'N1', N1
print >>h, 'MaxLen1', MaxLen1
print >>h, 'nHaps', 'NHaps', 'nSp', 'single_nucleotide_pDNCs', 'pDNC_recovered', 
print 'n', 'N', 'nSp','qPVS', 'TPVS','nPureDNCs', 'Retrieved_pDNCs', 'n_EDV_pDNCs'
for n in range(2, len(UHaplotypes[qCLADE])+1):#This block of code samples different number of species from all clades, this number is proportional to the number of species in a clade.
    N = 10
    while N > 0:
        dataset = open('D:\qHaps_pDNC_resampling\\dataset_'+qCLADE+'_'+gene+str(n)+str(N)+'.txt', 'w')
        raw_records = []
        for Clade in Clades1:
            m = int(n*(float(len(UHaplotypes[Clade])+0.1)) / len(UHaplotypes[qCLADE])) + 1
            if len(UHaplotypes[Clade]) > m:
                CHaps = random.sample(UHaplotypes[Clade], m)
            else:
                CHaps = UHaplotypes[Clade]
            for thing in CHaps:
                raw_records.append([Clade, Clade, thing])
                print >>dataset, Clade, Clade, thing
        dataset.close()
        
        Clades, clade_sorted_seqs, shared_positions = Step1(raw_records)#Functions in 1.1.
        SubSampleHaps = {}
        for smallclade in Clades:
            SubSampleHaps[smallclade] = Hcount(clade_sorted_seqs[smallclade])
            #float(DefineVarPos(SubSampleHaps[qCLADE]))/len(raw_records[0][2])
        print n, N, len(set([i[1] for i in raw_records])), float(len(DefineVarPos(SubSampleHaps[qCLADE])))/len(raw_records[0][2]), float(len(DefineVarPos([i[2] for i in raw_records])))/len(raw_records[0][2]),#len(raw_records),
        print >>h, n, N, len(set([i[1] for i in raw_records])), float(len(DefineVarPos(SubSampleHaps[qCLADE])))/len(raw_records[0][2]), float(len(DefineVarPos([i[2] for i in raw_records])))/len(raw_records[0][2]),
        x,y,z,pures = C_VP_PP(clade_sorted_seqs, qCLADE, shared_positions, Cutoff)#Functions in 1.2
        
        posrates = PosRatesDictNucl(y, SubSampleHaps[qCLADE], FCP, SCP, TCP, TrTable)
        zerorates = [key for key in posrates if len(posrates[key]) == 0]

        newy = {key:y[key] for key in y if not key in pures}
        ND_combinations = [[item] for item in pures]
        print len(pures),
        print >>h, len(pures),
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
#        print >>h, n, N, len(ND_combinations), len(Allpos)
        print len(ND_combinations),
        print >>h, len(ND_combinations),
        NValid = 0
        Validity_array = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
        for comb in ND_combinations:
            OnlyConserved = True
            for pos in comb:
                if pos in shared_positions1[qCLADE]:
                    continue
                else:
                    OnlyConserved = False
                    break
            Valid = False
            if OnlyConserved == True and ConditionD(comb, x1, y1) == True:
                Valid = True
                NValid += 1
            
            zeroesin = True
            if len([pos for pos in comb if pos in zerorates]) == 0:
                zeroesin = False
            
            if Valid == True and zeroesin == False:
                Validity_array[len(comb)-1][0] += 1
            elif Valid == False and zeroesin == False:
                Validity_array[len(comb)-1][1] += 1
            elif Valid == True and zeroesin == True:
                Validity_array[len(comb)-1][2] += 1
            else:
                Validity_array[len(comb)-1][3] += 1
        print NValid, str(Validity_array).replace('[[', '|').replace('], [', '|').replace(']]', '')
        print >>h, NValid, str(Validity_array).replace('[[', '|').replace('], [', '|').replace(']]', '')
        N-=1
h.close()
f.close()
