def euleriancycle(graph):
    stack = []
    circuit = []
    nodes = list(graph.keys())
    initial = nodes[0]
    current = nodes[0]
    while((len(graph[current])>0) | (len(stack)>0)):
        if(len(graph[current])<=0):
            circuit.append(current)
            current = stack.pop()
        else:
            stack.append(current)
            neighbour = graph[current][-1]
            graph[current].pop()
            current = neighbour
    circuit.append(initial)
    return circuit

graph = {0:[3],1:[0],2:[1,6],3:[2],4:[2],5:[4],6:[5,8],7:[9],8:[7],9:[6]}
path = euleriancycle(graph)
rev_path = path[::-1]
rp = []
for i in rev_path: #NOT req if graph is as str
    rp.append(str(i))
right_path = '->'.join(rp)


def calc_indeg(graph,node):
    indeg = 0
    for val in graph.values():
        #print(val)
        if ((len(val)>0)):
            for v in val:
                if(node == v):
                    indeg += 1
    return indeg


def calc_outdeg(graph,node):
    outdeg = 0
    #print(node)
    if node in graph.keys():
        outdeg = len(graph[node])
    return outdeg

def calc_inout(graph):
    all_nodes = []
    for key,val in graph.items():
        all_nodes.append(key)
        initial_node = key
        if(len(val)>1):
            for l in val:
                all_nodes.append(l)
        else:
            all_nodes.append(val)
    all_nodes = set(all_nodes)
    dict = {}
    for node in all_nodes:
        inout = []
        inout.append(calc_indeg(graph,node))
        inout.append(calc_outdeg(graph,node))
        dict[node] = inout
        if(inout[0]>inout[1]):
            initial_node = node
    return(initial_node)



def eulerianpath(graph):
    stack = []
    circuit = []
    all_nodes = []
    initial_node = list(graph.keys())[0]
    for key,val in graph.items():
        all_nodes.append(key)
        if(len(val)>0):
            for l in val:
                all_nodes.append(l)
    all_nodes = list(set(all_nodes))
    for n in all_nodes:
        if n not in graph.keys():
            graph[n]=[]
        ind = calc_indeg(graph,n)
        oud = calc_outdeg(graph,n)
        #print(f'node {n}: in: {ind}, out: {oud}')
        if(ind<oud):
            initial_node = n
    #print(initial_node)
    nodes = list(graph.keys())
    current = initial_node #calc_inout(graph)
    #stack.append(current)
    while((len(graph[current])>0) | (len(stack)>0)):
        #if (current not in graph.keys()):
        #    circuit.append(current)
        #    current = stack.pop()
        if(len(graph[current])<=0):
            circuit.append(current)
            current = stack.pop()
        else:
            stack.append(current)
            #print(stack)
            #for neigh in graph[current]:
            #    if (len(graph[neigh])>0):
            #        neighbour = neigh
            #        break
            #if(len(graph))
            neighbour = graph[current][0]
            graph[current].remove(neighbour)
            current = neighbour
    circuit.append(current)
    #print(stack)
    return circuit


{'0': ['1'], '1': ['2'], '2': ['3']}
graph = {0:[2],1:[3],2:[1],3:[0,4],6:[3,7],7:[8],8:[9],9:[6]}

with open('EulerianPath/inputs/test5.txt', 'r') as file:
with open('dataset_203_6.txt', 'r') as file:
    graph = dict((line.strip().split(' -> ') for line in file))
    for key in graph:
        graph[key] = graph[key].split(',')


path = eulerianpath(graph)
rev_path = path[::-1]
rp = []
for i in rev_path: #NOT req if graph is as str
    rp.append(str(i))


right_path = '->'.join(rp)
right_path


with open('dataset_203_6.txt', 'r') as file:
    graph = dict((line.strip().split(' -> ') for line in file))
    for key in graph:
        graph[key] = graph[key].split(',')

def stringreconstruction(patterns):
    db = debruijn(patterns)
    path = eulerianpath(db)
    rev_path = path[::-1]
    text = str(rev_path[0])
    for i in range(1, len(rev_path)):
        text = text+rev_path[i][-1]
    return (text)

GGCTTACCA
GGCGCTCTTTTATACACCCCA

pat = ['CTTA','ACCA','TACC','GGCT','GCTT','TTAC']
k = 4
stringreconstruction(k,pat)

with open('dataset_203_7.txt',"r") as file:
    text = file.read()
    klst = text.split('\n')
    genome = stringreconstruction(int(klst[0]),klst[1:-1])


def universal_str(k):
    kmers = ["".join(i) for i in itertools.product('01', repeat=k)]
    db = debruijn(kmers)
    path = euleriancycle(db)
    rev_path = path[::-1]
    text = str(rev_path[0])
    for i in range(1, len(rev_path)):
        text = text+rev_path[i][-1]
    return (text[0:-k+1])


g = 'TAATGCCATGGGATGTT'
g = 'TAATGCCATGGGATGTT'
ks = Compositionk(3,g)
l = []
for i in range(len(ks)-5):
    s = []
    s.append(ks[i])
    s.append(ks[i+5])
    l.append(s)

l = []
for i in range(len(ks)-5):
    s = "("
    s += ks[i]
    s += '|'
    s += ks[i+5]
    s += ')'
    l.append(s)

gappedpatterns = [['AG','AG'],['GC','GC'],['CA','CT'],['AG','TG'],['GC','GC'],['CT','CT'],['TG','TG'],['GC','GC'],['CT','CA']]
def stringspellecbygappedpatterns(gappedpatterns, k=3, d=3):
    fp = ''
    sp = ''
    for i in gappedpatterns:
        fp += i[0][0]
        sp += i[1][0]
    fp +=gappedpatterns[-1][0][1:]
    sp +=gappedpatterns[-1][1][1:]
    #print(fp)
    #print(sp)
    prefixstring = stringspellecbygappedpatterns(fp, k, d=3)


gp = [['GACC','GCGC'],['ACCG','CGCC'],['CCGA','GCCG'],['CGAG','CCGG'],['GAGC','CGGA']]

def conststr(kmers):
    fp = ''
    for i in kmers:
        fp +=i[0]
    fp += kmers[-1][1:]
    return fp

def pathgraphj(joined_kmers):
    nodes = joined_kmers
    lst = []
    current = nodes[0]
    lst.append(current)
    nodes.remove(current)
    while(len(nodes)>0):
        for j in nodes:
            if((j[0][0] == current[0][-1]) & (j[1][0] == current[1][-1])):
                print(lst)
                lst.append(j)
                nodes.remove(j)
                current = j
                #print(nodes)
    return lst





jp = [['TA', 'GC'], ['AT', 'GA'], ['TG', 'AT'], ['GG', 'TG'], ['GG', 'GT'], ['AA', 'CC'], ['AT', 'CA'], ['TG', 'AT'], ['GC', 'TG'], ['CC', 'GG'], ['CA', 'GG']]
jp2 = [['AA', 'CC'], ['TG', 'AT'], ['GG', 'TG'], ['GG', 'GT'], ['GA', 'TT'], ['AT', 'CA'], ['TG', 'AT'], ['GC', 'TG'], ['CC', 'GG'], ['CA', 'GG'], ['AT', 'GA']]
def stringspellecbygappedpatterns(gappedpatterns, k, d):
    fp = []
    sp = []
    for i in gappedpatterns:
        fp.append(i[0])
        sp.append(i[1])
        #fp.append([i[0][:k-1],i[1][:k-1]])
        #sp.append([i[0][1:],i[1][1:]])
    #print(fp)
    #print("----------------")
    #print(sp)
    prefixstring =conststr(fp)#stringreconstruction(k,fp) # conststr(fp) #
    #print(prefixstring)
    #print("---------------------------------------------")
    suffixstring =conststr(sp)#stringreconstruction(k,sp) # conststr(sp)#
    #print(suffixstring)
    for i in range(k+d, len(prefixstring)):
        if(prefixstring[i] != suffixstring[i-k-d]):
            return ('No string spelled by gapped pattern')
    return(prefixstring+suffixstring[-(d+k):])
gp = [['GAGA','TTGA'],['TCGT','GATG'],['CGTG','ATGT'],['TGGT','TGAG'],['GTGA','TGTT'],['GTGG','GTGA'],['TGAG','GTTG'],['GGTC','GAGA'],['GTCG','AGAT']]


def pairedreads(gp):
    dict = {}
    for i in gp:
        #print(i)
        key = i[0][:-1]+i[1][:-1]
        value = i[0][1:]+i[1][1:]
        if(dict.get(key,0)==0):
            dict[key] = [value]
        else:
            dict[key].append(value)
    return dict





def qpairread(gp,k,d):
    dd = pairedreads(gp)
    path = eulerianpath(dd)
    rev_path = path[::-1]
    l = []
    for i in range(len(rev_path)-1):
        j = i+1
        s1 = rev_path[i][:k-1] + rev_path[j][k-2]
        s2 = rev_path[i][k-1:] + rev_path[j][-1]
        l.append([s1,s2])
    text = stringspellecbygappedpatterns(l, k, d)
    return text


#print(l)
print('===========================================================================')
stringspellecbygappedpatterns(l, k, d)




gp = [['TAA','GCC'],['ATG','GAT'],['TGG','ATG'],['GGG','TGT'],['GGA','GTT'],['AAT','CCA'],['ATG','CAT'],['TGC','ATG'],['GCC','TGG'],['CCA','GGG'],['CAT','GGA']]
k=3
d=1

prefixes, suffixes = [stringreconstruction(patterns) for p in zip(*[(node[:k-1], node[k-1:]) for node in gp])]

with open('PairedStringReconstruction/inputs/test3.txt', 'r') as file:

with open('dataset_204_16.txt', 'r') as file:
with open('reads.txt', 'r') as file:
    text = file.read()
    lst = text.split('\n')
    #k,d = lst[0].split()
    #k = int(k)
    #d = int(d)
    gp = list(line.strip().split('|') for line in lst[:-1])

k= 120
d = 1000
t = qpairread(gp,k,d)
t

dd = pairedreads(gp)
path = eulerianpath(dd)
rev_path = path[::-1]
l = []
for i in range(len(rev_path)-1):
    j = i+1
    s1 = rev_path[i][:k-1] + rev_path[j][k-2]
    s2 = rev_path[i][k-1:] + rev_path[j][-1]
    l.append([s1,s2])


#print(l)
print('===========================================================================')
stringspellecbygappedpatterns(l, k, d)



text = stringspellecbygappedpatterns(gp, int(k), int(d))



g = {1:[2], 2:[3], 3:[4,5], 6:[7], 7:[6]}

def maximalnonbranchingpaths(graph):
    paths = []
    joined = []
    for v in graph.keys():
        if not ((calc_indeg(graph,v) == 1) & (calc_outdeg(graph,v) == 1)):
            print('in loop 1')
            if(calc_outdeg(graph,v) > 0):
                nonbranchingpath = [v]
                for w in graph[v]:
                    nonbranchingpath = [v]
                    nonbranchingpath.append(w)
                    joined.append(v)
                    joined.append(w)
                    while ((calc_indeg(graph,w) == 1) & (calc_outdeg(graph,w) == 1)):
                        nonbranchingpath.append(graph[w][0])
                        joined.append(w)
                        joined.append(graph[w][0])
                        w = graph[w][0]
                        print(nonbranchingpath)
                    paths.append(nonbranchingpath)
    #path = eulerianpath(graph)
    for v in graph.keys():
        if ((calc_indeg(graph,v) == 1) & (calc_outdeg(graph,v) == 1)):
            print('in loop 2')
            if(v not in joined):
                nb = [v]
                for w in graph[v]:
                    nb.append(w)
                    joined.append(v)
                    joined.append(w)
                    while ((calc_indeg(graph,w) == 1) & (calc_outdeg(graph,w) == 1)& (w !=v)):
                        nb.append(graph[w][0])
                        joined.append(w)
                        w = graph[w][0]
                    paths.append(nb)
    return (paths)

    #    nb = []
    #    if((len(v)>0) & (calc_indeg(graph,k) == 1) & (calc_outdeg(graph,k) == 1)):
    #        nb.append([k,v])
    #        while ((calc_indeg(graph,v) == 1) & (calc_outdeg(graph,v) == 1)):
    #            nb.append([v,graph[v]])
    #            v = graph(v)
    #        paths.append(nb)
    #return paths

with open('MaximalNonBranchingPaths/inputs/test5.txt', 'r') as file:
with open('dataset_6207_2.txt', 'r') as file:
    graph = dict((line.strip().split(' -> ') for line in file))
    for key in graph:
        graph[key] = graph[key].split(',')

cycles = maximalnonbranchingpaths(graph)

for c in cycles:
    right_path = '->'.join(c)
    print(right_path)

kmers = ['ATG','ATG','TGT','TGG','CAT','GGA','GAT','AGA']

with open('dataset_205_5.txt',"r") as file:
    text = file.read()
    klst = text.split('\n')


['AAAT','AATG','ACCC','ACGC','ATAC','ATCA','ATGC','CAAA','CACC','CATA','CATC','CCAG','CCCA','CGCT','CTCA','GCAT','GCTC','TACG','TCAC','TCAT','TGCA']

gp = [['ACC','ATA'],['ACT','ATT'],['ATA','TGA'],['ATT','TGA'],['CAC','GAT'],['CCG','TAC'],['CGA','ACT'],['CTG','AGC'],['CTG','TTC'],['GAA','CTT'],['GAT','CTG'],['GAT','CTG'],['TAC','GAT'],['TCT','AAG'],['TGA','GCT'],['TGA','TCT'],['TTC','GAA']]

d = debruijn(klst[:-1])
cy = maximalnonbranchingpaths(d)
l = []
for i in cy:
    s = i[0]
    for code in i[1:]:
        s += code[-1]
    l.append(s)


l.sort()
print(*l)
