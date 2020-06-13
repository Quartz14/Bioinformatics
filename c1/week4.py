import random
def randomizedmotifsearch(dna,k,t):
    motifs = []
    for line in dna:
        kmers = [line[i:i+k] for i in range(len(line)-k+1)]
        index = random.randint(0,len(kmers)-1)
        motifs.append(kmers[index])
    best_motifs = motifs[:]
    #print(best_motifs)
    while(True):
        profile = laplace_profile(motifs,k,t,1)
        #print(profile)
        motifs = []
        for i in range(t):
             m = profile_most_prob(dna[i],k,profile)
             motifs.append(m)
        #motifs = profile_to_motif(dna,k,profile)
        if(score_diff(motifs,k,t)< score_diff(best_motifs,k,t)):
            best_motifs = motifs[:]
        else:
            return(best_motifs)
motifs = ['TAAC','GTCT','CCGG','ACTA','AGGT']
k=4
t=5
dna = ['TTACCTTAAC','GATGTCTGTC','CCGGCGTTAG','CACTAACGAG','CGTCAGAGGT']
profile = laplace_profile(motifs,k,t,1)
for i in range(t):
     m = profile_most_prob(dna[i],k,profile)
     print(m)


dna = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA','GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG','TAGTACCGAGACCGAAAGAAGTATACAGGCGT','TAGATCAAGTTTCAGGTGCACGTCGGTGAACC','AATCCACCAGCTCCACGTGCAATGTTGGCCTA']


def randomizedmotifsearch_loop(dna,k,t):
    best_motifs =  randomizedmotifsearch(dna,k,t)
    for i in range(1000):
        motifs =  randomizedmotifsearch(dna,k,t)
        if(score_diff(motifs,k,t)< score_diff(best_motifs,k,t)):
            print(score_diff(motifs,k,t))
            best_motifs = motifs[:]
    return(best_motifs)

print(best_motifs)



with open("dataset_161_5.txt",'r') as file:
    lines = file.read().split('\n')
    k,t = lines[0].split()
    #print(k,d)
    dna = lines[1:-1]
    #print(dna)
    best_motifs2 =  randomizedmotifsearch(dna,int(k),int(t))
    for i in range(1000):
        motifs =  randomizedmotifsearch(dna,int(k),int(t))
        if(score_diff(motifs,int(k),int(t))< score_diff(best_motifs,int(k),int(t))):
            #print(score_diff(motifs,int(k),int(t)))
            best_motifs2 = motifs[:]
    print(*best_motifs2,sep='\n')

def gibbssampler(dna,k,t,n):
    motifs = []
    for line in dna:
        kmers = [line[i:i+k] for i in range(len(line)-k+1)]
        index = random.randint(0,len(kmers)-1)
        motifs.append(kmers[index])
    best_motifs = motifs[:]
    for j in range(n):
        i = random.randint(0,t-1)
        sub_motif = []
        for id,m in enumerate(motifs):
            if(id!=i):
                sub_motif.append(m)
        profile = laplace_profile(sub_motif,k,t-1,1)
        motifs[i] = profile_random_prob(dna[i],k,profile)
        #print(motifs)
        if(score_diff(motifs,k,t)< score_diff(best_motifs,k,t)):
            #print(score_diff(motifs,k,t))
            best_motifs = motifs[:]
    return(best_motifs)

dna = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA','GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG','TAGTACCGAGACCGAAAGAAGTATACAGGCGT','TAGATCAAGTTTCAGGTGCACGTCGGTGAACC','AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
k = 8
t = 5
n = 100
best_motifs =  gibbssampler(dna,k,t,n)
for i in range(20):
    motifs =  gibbssampler(dna,k,t,n)
    if(score_diff(motifs,int(k),int(t))< score_diff(best_motifs,int(k),int(t))):
        best_motifs = motifs[:]
print(*best_motifs,sep='\n')

with open("dataset_163_4.txt",'r') as file:
    lines = file.read().split('\n')
    k,t,n = lines[0].split()
    #print(k,t,n)
    dna = lines[1:-1]
    #print(dna)
    k = int(k)
    t= int(t)
    n = int(n)
    best_motifs =  gibbssampler(dna,k,t,n)
    for i in range(20):
        motifs =  gibbssampler(dna,k,t,n)
        if(score_diff(motifs,int(k),int(t))< score_diff(best_motifs,int(k),int(t))):
            #print(score_diff(motifs,int(k),int(t)))
            best_motifs = motifs[:]
    print(*best_motifs,sep='\n')


with open("DosR.txt",'r') as file:
    lines = file.read().split('\n')
    #k,t,n = lines[0].split()
    #print(k,t,n)
    dna = lines[:-1]
    #print(dna)
    #print(lines)
    kss = [8,9,10,11,12]
    t=10
    for k in kss:
        best_motifs2 =  randomizedmotifsearch(dna,int(k),int(t))
        for i in range(1000):
            motifs =  randomizedmotifsearch(dna,int(k),int(t))
            if(score_diff(motifs,int(k),int(t))< score_diff(best_motifs,int(k),int(t))):
                #print(score_diff(motifs,int(k),int(t)))
                best_motifs2 = motifs[:]
        print(*best_motifs2,sep='\n')
        print('-----------------------------------')




def laplace_profile(dna,k,t,pseudocount=1):
    a = []
    c = []
    tt = []
    g = []
    for i in range(len(dna[0])):
        a.append(pseudocount)
        c.append(pseudocount)
        tt.append(pseudocount)
        g.append(pseudocount)
    profile = {'A':a,'C':c,'G':g,'T':tt}
    for kmer in dna:
        for j,char in enumerate(kmer):
            profile[char][j]+=1
    for val in profile.values():
        for i,num in enumerate(val):
            val[i] = val[i]/(t+(4*pseudocount))
    return profile

def profile_to_motif(dna,k,profile):
    motifs = []
    for line in dna:
        kmers = [line[i:i+k] for i in range(len(line)-k+1)]
        probabilities = []
        for pattern in kmers:
            prob = 1
            letter=0
            for i,char in enumerate(pattern):
                prob *= profile[char][i]
            probabilities.append(prob)
        maxprob = max(probabilities)
        index = probabilities.index(max(probabilities))
        motifs.append(kmers[index])
    return motifs

def score_diff(Motif,k,t): #as dict
    concenses=''
    profile = laplace_profile(Motif,k,t,1)
    for i in range(k):
        max_val = 0
        for key in profile.keys():
            val = profile[key][i]
            if(val>max_val):
                max_val = val
                char = key
        #print(char)
        concenses += char
    #print(concenses)
    diff = 0
    for mer in Motif:
        diff += hamming_dist(concenses, mer)
    return(diff)

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

import numpy as np
def profile_random_prob(dna,k,profile):   #profile is dict, dna is ith strand
    kmers = [dna[i:i+k] for i in range(len(dna)-k+1)]
    probabilities = []
    for pattern in kmers:
        prob = 1
        letter=0
        for i,char in enumerate(pattern):
            prob *= profile[char][i]
        probabilities.append(prob)
    #randomprob = random.choice(probabilities)
    #index = probabilities.index(randomprob)
    psum = 0
    for p in probabilities:
        psum +=p
    updated_prob = []
    for prob in probabilities:
        updated_prob.append(prob/psum)
    index = random.choice(range(len(probabilities)))
    random_k = np.random.choice(kmers, 1,True,updated_prob)
    #print(random_k[0])
    return random_k[0]
def hamming_dist(text1, text2):
    n1 = len(text1)
    n2 = len(text2)
    dist = 0
    if(n1==n2):
        for i in range(n1):
            if(text1[i] != text2[i]):
                dist +=1
    return dist



motifs = ['GTC','CCC','ATA','GCT']
k=3
t=4
profile = laplace_profile(motifs,k,t,1)
#print(profile)
motifs = []
for i in range(t):
     m = profile_most_prob(dna[i],k,profile)
     motifs.append(m)

motifs = ['TGA','GTT','GAA','TGT']
k=3
t=4
dna = ['TGACGTTC','TAAGAGTT','GGACGAAA','CTGTTCGC']
while(True):
    profile = laplace_profile(motifs,k,t,1)
    #print(profile)
    motifs = []
    for i in range(t):
         m = profile_most_prob(dna[i],k,profile)
         motifs.append(m)
    #motifs = profile_to_motif(dna,k,profile)
    print(motifs)
    if(score_diff(motifs,k,t)< score_diff(best_motifs,k,t)):
        best_motifs = motifs[:]
    else:
        print(*best_motifs)




TTCGTGACCGACGTCCCCAG
TTGGGGACTTCCGGCCCTAA
GCCGGGACTTCAGGCCCTAT
CATGGGACTTTCGGCCCTGT
GAGGGGACTTTTGGCCACCG
CCAGGGACCTAATTCCATAT
TTGAGGACCTTCGGCCCCAC
CTGGGGACCGAAGTCCCCGG
TTAGGGACCATCGCCTCCTG
TGGATGACTTACGGCCCTGA
TTGGGGACTAAAGCCTCATG
TCGGGGACTTCTGTCCCTAG
TTGGGGACCATTGACCCTGT
TTGAGGACCTAAGCCCGTTG
CACGGGTCAAACGACCCTAG
GGCGGGACGTAAGTCCCTAA
GAAGTGACGAAAGACCCCAG
CGGAGGACCTTTGGCCCTGC
GTGGGGACCAACGCCCCTGG
