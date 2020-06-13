s = "CAATCCAAC"

def Compositionk(k,text):
    kmers = []
    for i in range(len(text)-k+1):
        kmers.append(text[i:i+k])
    return kmers

ks = Compositionk(5,s)
print(*ks)

with open('dataset_197_3.txt',"r") as file:
    text = file.read()
    k, s,_ = text.split('\n')
    kms = Compositionk(int(k),s)

kms.sort()
for i in kms:
    print(i)

txt = '\n'.join(kms)
f = open('kms_sol.txt','w')
f.write(txt)
f.close()

if __name__ == "__main__":
    Input = sys.stdin.readlines()
    patternList = [pattern.strip() for pattern in Input]
    ans = Pathtogenome(patternList)
    print(ans)

def Pathtogenome(kmers_lst):
    l = len(kmers_lst)
    genome=''
    if (l<=0):
        return genome
    for i in range(l-1):
        genome = genome+kmers_lst[i][0]
    genome = genome+kmers_lst[-1]
    return genome

klist = ['ACCGA','CCGAA','CGAAG','GAAGC','AAGCT']
Pathtogenome(klist)

with open('dataset_198_3.txt',"r") as file:
    text = file.read()
    klst = text.split('\n')
    genome = Pathtogenome(klst[:-1])




def overlap(kmers):
    dict = {}
    for i in kmers:
        dict[i] = []
    for i in range(len(kmers)):
        for j in range(len(kmers)):
            if((kmers[i] != kmers[j])&(kmers[i][1:]==kmers[j][:-1])):
                dict[kmers[i]].append(kmers[j])
    return {k: v for k, v in dict.items() if len(v)>0}

kmers = ['ATGCG','GCATG','CATGC','AGGCA','GGCAT','GGCAC']
graph = overlap(kmers)

for k,v in graph.items():
    vlist = ", ".join(v)
    print(f'{k} -> {vlist}')

with open('dataset_198_10.txt',"r") as file:
    text = file.read()
    klst = text.split('\n')
    graph = overlap(klst[:-1])

def pathgraph(k, genome):
    dict = {}
    kmers_lst = Compositionk(k, genome)
    #print(genome[0:k-1],genome[1:])
    #dict[genome[0:k-1]] = []
    for i in kmers_lst:
        dict[i[:-1]] = []
    dict[genome[-(k-1):]]= []
    for i in range(len(kmers_lst)-1):
        if(genome[0:k-1] == kmers_lst[i][:-1]):
            dict[genome[0:k-1]].append(kmers_lst[i][1:])
        #for j in range(i+1,len(kmers_lst)):
        j = i+1
        if((kmers_lst[i] != kmers_lst[j]) & (kmers_lst[i][1:] == kmers_lst[j][:-1])):
            dict[kmers_lst[i][1:]].append(kmers_lst[j][1:])
    return {k: v for k, v in dict.items() if len(v)>0}

d = pathgraph(3,'TAATGCCATGGGATGTT')
d = pathgraph(4,'AAGATTCTCTAAGA')
for k,v in d.items():
    v.sort()
d

with open('dataset_199_6.txt',"r") as file:
    text = file.read()
    k, s,_ = text.split('\n')
    d = pathgraph(int(k),s)

for k,v in d.items():
    v.sort()

keys = list(d.keys())
keys.sort()
for k in keys:
    vlist = ", ".join(d[k])
    print(f'{k} -> {vlist}')

def debruijn(patterns):
    dict = {}
    unique_ps = []
    for kmer in patterns:
        unique_ps.append(kmer[:-1])
        unique_ps.append(kmer[1:])
    unique_ps = set(unique_ps)
    for p in unique_ps:
        dict[p] = []
    for kmer in patterns:
        dict[kmer[:-1]].append(kmer[1:])
    return {k: v for k, v in dict.items() if len(v)>0}

l = ['GCGA','CAAG','AAGA','GCCG','ACAA','AGTA','TAGG','AGTA','ACGT','AGCC','TTCG','AGTT','AGTA','CGTA','GCGC','GCGA','GGTC','GCAT','AAGC','TAGA','ACAG','TAGA','TCCT','CCCC','GCGC','ATCC','AGTA','AAGA','GCGA','CGTA']

['GAGG','CAGG','GGGG','GGGA','CAGG','AGGG','GGAG']
d = debruijn(['GAGG','CAGG','GGGG','GGGA','CAGG','AGGG','GGAG'])

with open('dataset_200_8.txt',"r") as file:
    text = file.read()
    klst = text.split('\n')

print_graph(['GCAAG','CAGCT','TGACG']) #TEST 1
GCAA -> CAAG
CAGC -> AGCT
TGAC -> GACG
print_graph(['AGGT','GGCT','AGGC']) #TEST 2
AGG -> GGT,GGC
GGC -> GCT
print_graph(['TTCT','GGCT','AAGT','GGCT','TTCT']) #TEST 3
TTC -> TCT,TCT
GGC -> GCT,GCT
AAG -> AGT
print_graph(['CA','CA','CA','CA','CC','CA']) #TEST 4
C -> A,A,A,A,C,A
lst = ['GAGG','CAGG','GGGG','GGGA','CAGG','AGGG','GGAG']
print_graph(lst)

print_graph(klst[:-1])
def print_graph(patterns):
    d = debruijn(patterns)
    #for k,v in d.items():
    #    v.sort()
    #keys = list(d.keys())
    #keys.sort()
    #for k in keys:
    for k,v in d.items():
        vlist = ",".join(v)
        print(f'{k} -> {vlist}')


d = debruijn(klst[:-1])

for k,v in d.items():
    v.sort()

keys = list(d.keys())
keys.sort()

for k in keys:
    vlist = ", ".join(d[k])
    print(f'{k} -> {vlist}')
