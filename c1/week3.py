def motifenumeration(dna, k, d):
    patterns = []
    kmers = []
    neighborhood = []
    motif_dict = {}
    #for string in dna:
    string = dna[0]
    no_of_lines = len(dna)
    for i in range(len(string)-k+1):
        kmers.append(string[i:i+k])
        neighborhood.append(Neighbors(kmers[i],d))
    #print(len(kmers))
    #print(len(neighborhood))
    #print(no_of_lines)
    for i in range(len(kmers)):
        motif_dict[kmers[i]] = neighborhood[i]
    #print(motif_dict)
    for pattern,pattern_ in motif_dict.items():
        for pat_ in pattern_:
            count = 0
            for string in dna:
                for i in range(len(string)-k+1):
                    if( hamming_dist(pat_, string[i:i+k])<=d):
                        count +=1
                        break
            if(count==no_of_lines):
                patterns.append(pat_)
    return(set(patterns))

with open("dataset_156_8.txt",'r') as file:
    lines = file.read().split('\n')
    k,d = lines[0].split()
    #print(k,d)
    dna = lines[1:-1]
    #print(dna)
    patterns = motifenumeration(dna, int(k),int(d))
    print(*patterns)


motifenumeration(['ATTTGGC','TGCCTTA','CGGTATC','GAAAATT'],3,1)

Motifs =["TCGGGGGTTTTT","CCGGTGACTTAC","ACGGGGATTTTC","TTGGGGACTTTT","AAGGGGACTTCC","TTGGGGACTTCC","TCGGGGATTCAT","TCGGGGATTCCT","TAGGGGAACTAC","TCGGGTATAACC"]

import math
def motif_sum(profile):
    final=0
    for col in profile:
        #print(col)
        sum = 0
        for val in col:
            if(val !=0):
                sum += val*math.log(val,2)
        #print(sum)
        final += sum
    return(final)

n = len(profile['A'])

def score(profile,k): #as dict
    final=0
    for i in range(k):
        sum = 0
        for key in profile.keys():
            val = profile[key][i]
            if(val !=0):
                sum += val*math.log(val,2)
        final +=sum
    return (-1*final)

def score_diff(Motif,k,t): #as dict
    concenses=''
    profile = profile_mat(Motif,k,t)
    for i in range(k):
        max_val = 0
        for key in profile.keys():
            val = profile[key][i]
            if(val>max_val):
                max_val = val
                char = key
        concenses += char
    #print(concenses)
    diff = 0
    for mer in Motif:
        diff += hamming_dist(concenses, mer)
    return(diff)

profile ={'A':[0.4 ,0.3 ,0.0 ,0.1 ,0.0 ,0.9],'C': [0.2 ,0.3 ,0.0 ,0.4 ,0.0 ,0.1],'G': [0.1 ,0.3 ,1.0 ,0.1 ,0.5 ,0.0],'T': [0.3 ,0.1 ,0.0 ,0.4 ,0.5 ,0.0]}
d = ['GGC','AAG','CAA','CAC','CAA']
k=3
t=5
d = ['GCG','AAG','AAG','ACG','CAA']
d = ['CGT','AAG','AAG','AAT','AAT']
d = ['CAG','CAG','CAA','CAA','CAA']
score_diff({'A': [0.4, 0.8, 0.2], 'C': [0.4, 0.0, 0.4], 'G': [0.2, 0.2, 0.4], 'T': [0.0, 0.0, 0.0]},3)
score({'A': [0.4, 0.8, 0.2], 'C': [0.4, 0.0, 0.4], 'G': [0.2, 0.2, 0.4], 'T': [0.0, 0.0, 0.0]},3)
for key,val in profile.items:
    word = profile[key]

#Do not enter 0 as its undefined
profile = [[0.2,0.1,0.7],[0.2,0.6,0.2],[1]
            ,[1],[0.9,0.1],[0.9,0.1],
            [0.9,0.1],[0.1,0.4,0.5],[0.1,0.1,0.8],
            [0.1,0.2,0.7],[0.3,0.4,0.3],[0.6,0.4]]



def distance_bwt_pattern_string(pattern,dna):
    k = len(pattern)
    distance=0
    for text in dna:
        text_kmer = []
        hammingdist = float('inf')
        for i in range(len(text)-k+1):
            text_kmer.append(text[i:i+k])
        for pattern_ in text_kmer:
            hd = hamming_dist(pattern,pattern_)
            if(hammingdist>hd):
                hammingdist = hd
        distance = distance+hammingdist
    return(distance)

distance_bwt_pattern_string('AAA',['TTACCTTAAC', 'GATATCTGTC', 'ACGGCGTTCG', 'CCCTAAAGAG', 'CGTCAGAGGT'])

with open("dataset_5164_1.txt",'r') as file:
    lines = file.read().split('\n')
    pattern = lines[0]
    #print(pattern)
    dna = lines[1].split()
    #print(dna)
    distance = distance_bwt_pattern_string(pattern,dna)
    print(distance)

def medianstring(dna,k):
    distance = float('inf')
    for i in range(4**k):
        pattern = numbertopattern(i,k)
        dbps = distance_bwt_pattern_string(pattern,dna)
        if (distance >dbps ):
            distance =  dbps
            Median = pattern
    #print(Median)
    return Median
dna = ['CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC','GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC','GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG']
k = 7
dna = ["AAATTGACGCAT", "GACGACCACGTT", "CGTCAGCGCCTG", "GCTGAGCACCGG", "AGTTCGGGACAG"]
medianstring(dna,3)
with open("dataset_158_9.txt",'r') as file:
    lines = file.read().split('\n')
    #print(lines)
    k = lines[0]
    #print(k)
    dna = lines[1:-1]
    #print(dna)
    median = medianstring(dna,int(k))
    print(median)

dna = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
profile = {
    'A': [0.2, 0.2, 0.3, 0.2, 0.3],
    'C': [0.4, 0.3, 0.1, 0.5, 0.1],
    'G': [0.3, 0.3, 0.5, 0.2, 0.4],
    'T': [0.1, 0.2, 0.1, 0.1, 0.2]}
k = 5
def profile_most_prob(dna,k,profile):   #profile is dict
    kmers = [dna[i:i+k] for i in range(len(dna)-k+1)]
    probabilities = []
    for pattern in kmers:
        prob = 1
        letter=0
        for i,char in enumerate(pattern):
            prob *= profile[char][i]
        probabilities.append(prob)
    maxprob = max(probabilities)
    index = probabilities.index(max(probabilities))
    return kmers[index]

profile_most_prob(dna,k,profile)
with open("dataset_159_3.txt",'r') as file:
    lines = file.read().split('\n')
    #print(lines)
    dna = lines[0]
    k = lines[1]
    #print(k)
    profile = lines[2:-1]
    #print(dna)
    clean = {}
    letters = ['A','C','G','T']
    for i,line in enumerate(profile):
        p = line.split()
        p = [float(i) for i in p]
        clean[letters[i]] = p
    #print(clean)
    most_prob = profile_most_prob(dna,int(k),clean)
    print(most_prob)

dna = ['GGCGTTCAGGCA','AAGAATCAGTCA','CAAGGAGTTCGC','CACGTCAATCAC','CAATAATATTCG']
k = 3
t = 5

def greedymotifsearch(dna,k,t):
    best_motifs = []
    for line in dna:
        best_motifs.append(line[:k])
    #print(best_motifs)
    kmers = [dna[0][i:i+k] for i in range(len(dna[0])-k+1)]
    #print(kmers)
    for motif in kmers:
        Motif = [motif]
        for i in range(1,t):
            profile = profile_mat(Motif,k,t)
            most_prob = profile_most_prob(dna[i],k,profile)
            Motif.append(most_prob)
        if(score(profile_mat(Motif,k,t),k) <score(profile_mat(best_motifs,k,t),k) ):
            best_motifs = Motif
    return(best_motifs)

def greedymotifsearchv2(dna,k,t):
    best_motifs = []
    for line in dna:
        best_motifs.append(line[:k])
    #print(best_motifs)
    kmers = [dna[0][i:i+k] for i in range(len(dna[0])-k+1)]
    #print(kmers)
    for motif in kmers:
        Motif = [motif]
        for i in range(1,t):
            profile = profile_mat(Motif,k,t)
            most_prob = profile_most_prob(dna[i],k,profile)
            Motif.append(most_prob)
        if(score_diff(Motif,k,t) <score_diff(best_motifs,k,t) ):
            best_motifs = Motif
    return(best_motifs)

with open("dataset_159_5.txt",'r') as file:
    lines = file.read().split('\n')
    k,d = lines[0].split()
    #print(k,d)
    dna = lines[1:-1]
    #print(dna)
    patterns2 = greedymotifsearchv2(dna, int(k),int(d))
    print(*patterns2)

def greedymotifsearchlaplace(dna,k,t):
    best_motifs = []
    for line in dna:
        best_motifs.append(line[:k])
    #print(best_motifs)
    kmers = [dna[0][i:i+k] for i in range(len(dna[0])-k+1)]
    #print(kmers)
    for motif in kmers:
        Motif = [motif]
        for i in range(1,t):
            profile = laplace_profile(Motif,k,t)
            most_prob = profile_most_prob(dna[i],k,profile)
            Motif.append(most_prob)
        if(score(laplace_profile(Motif,k,t),k) <score(laplace_profile(best_motifs,k,t),k) ):
            best_motifs = Motif
    return(best_motifs)
with open("dataset_160_9.txt",'r') as file:
    lines = file.read().split('\n')
    k,d = lines[0].split()
    #print(k,d)
    dna = lines[1:-1]
    #print(dna)
    patterns = greedymotifsearchlaplace(dna, int(k),int(d))
    print(*patterns)

def laplace_profile(dna,k,t):
    a = []
    c = []
    tt = []
    g = []
    for i in range(k):
        a.append(1)
        c.append(1)
        tt.append(1)
        g.append(1)
    profile = {'A':a,'C':c,'G':tt,'T':g}
    for kmer in dna:
        for j,char in enumerate(kmer):
            profile[char][j]+=1
    for val in profile.values():
        for i,num in enumerate(val):
            val[i] = val[i]/(t+4)
    return profile


def profile_mat(dna,k,t):
    #['GAT','TAG','ACT']
    a = []
    c = []
    tt = []
    g = []
    for i in range(k):
        a.append(0)
        c.append(0)
        tt.append(0)
        g.append(0)
    profile = {'A':a,'C':c,'G':tt,'T':g}
    #print(profile)
    for kmer in dna:
        for j,char in enumerate(kmer):
            profile[char][j]+=1
            #print(char,profile[char])
    for val in profile.values():
        #print(val)
        for i,num in enumerate(val):
            #print(type(val[i]))
            #print(type(t))
            val[i] = val[i]/t
            #print(i,num)
            #print(val[i]/10)
    return profile

dna =["TCGGGGGTTTTT","CCGGTGACTTAC","ACGGGGATTTTC","TTGGGGACTTTT","ATGGGGACTTCC","TCGGGGACTTCC","TCGGGGATTCAT","TAGGGGATTCCT","TAGGGGAACTAC","TCGGGTATAACC"]
profile_mat(dna,12,10)




def Neighbors(Pattern, d):
    if(d==0):
        return([Pattern])
    if(len(Pattern)==1):
        return(['A','C','G','T'])
        #return(pattern)
    neighborhood = []
    suffixneighbors = Neighbors(Pattern[1:],d)
    #print(suffixneighbors)
    for text in suffixneighbors:
        if(hamming_dist(Pattern[1:],text)<d):
            for nucleotide in ['A','C','G','T']:
                neighborhood.append(nucleotide+text)
        else:
            neighborhood.append(Pattern[0]+text)
    return(neighborhood)

def hamming_dist(text1, text2):
    n1 = len(text1)
    n2 = len(text2)
    dist = 0
    if(n1==n2):
        for i in range(n1):
            if(text1[i] != text2[i]):
                dist +=1
    return dist
def numbertopattern(index, k):
    index = int(index)
    k = int(k)
    if (k == 1):
        return (numbertosymbol(index))
    prefixIndex = index//4
    r = index%4
    symbol = numbertosymbol(r)
    prefixpattern = numbertopattern(prefixIndex,k-1)
    return (str(prefixpattern)+str(symbol))

def numbertosymbol(number):
    if (number==0):
        return('A')
    elif(number==1):
        return('C')
    elif(number==2):
        return('G')
    elif(number==3):
        return('T')
