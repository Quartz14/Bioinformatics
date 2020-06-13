def skew(text):
    l = [0]
    res=[]
    for char in text:
        if(char == "C"):
            l.append(l[-1]-1)
        elif(char == "G"):
            l.append(l[-1]+1)
        else:
            l.append(l[-1])
    min_val = min(l)
    for i,val in enumerate(l):
        if(val==min_val):
            res.append(i)
    return(l,res)

l = skew("CATGGGCATCGGCCATACGCC")
print(*l)

with open('dataset_7_6.txt',"r") as file:
    text = file.read()
    l = skew(text)
    print(*l)

def hamming_dist(text1, text2):
    n1 = len(text1)
    n2 = len(text2)
    dist = 0
    if(n1==n2):
        for i in range(n1):
            if(text1[i] != text2[i]):
                dist +=1
    return dist

d = hamming_dist('GGGCCGTTGGT','GGACCGTTGAC')
print(d)

with open('dataset_9_3.txt',"r") as file:
    text1,text2,_= file.read().split('\n')
    d = hamming_dist(text1,text2)
    print(d)

def approx_match(pattern,text,d):
    l = []
    for i in range(len(text)-len(pattern)+1):
        subtext = text[i:i+len(pattern)]
        if(hamming_dist(subtext,pattern)<=d):
            l.append(i)
    return(l)

d = approx_match('ATTCTGGA','CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT',3)
print(*d)

with open('dataset_9_4.txt',"r") as file:
    pattern,text,d,_= file.read().split('\n')
    l = approx_match(pattern,text,int(d))
    print(*l)

d = approx_match('AAAAA','AACAAGCTGATAAACATTTAAAGAG',2)
print(d)
print(len(d))
d = approx_match('GAGG','TTTAGAGCCTTCAGAGG',2)
print(len(d))
with open('dataset_9_6.txt',"r") as file:
    pattern,text,d,_= file.read().split('\n')
    l = approx_match(pattern,text,int(d))
    print(len(l))

def symboltonumber(symbol):
    if (symbol=='A'):
        return(0)
    elif(symbol=='C'):
        return(1)
    elif(symbol=='G'):
        return(2)
    elif(symbol=='T'):
        return(3)

def patterntonumber(pattern):
    if len(pattern) == 0:
        return (0)
    symbol = pattern[-1]
    prefix = pattern[:-1]
    return (4*patterntonumber(prefix) + symboltonumber(symbol))

def numbertosymbol(number):
    if (number==0):
        return('A')
    elif(number==1):
        return('C')
    elif(number==2):
        return('G')
    elif(number==3):
        return('T')

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

def Neighbors(Pattern, d):
    if(d==0):
        return(list(Pattern))
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

def computingFrequencywithmissmatches(text,k,d):
    frequencyArray = []
    result = []
    neighborhood = []
    for i in range(4**k):
        frequencyArray.append(0)
    for i in range(len(text)-k+1):
        #pattern = text[i:i+k]
        neighborhood.append(Neighbors(text[i:i+k],d))
        #neighborhood = Neighbors(pattern,d)
    neighborhood = [j for l1 in neighborhood for j in l1]
    for str in neighborhood:
        j = patterntonumber(str)
        frequencyArray[j] +=1
    maxcount = max(frequencyArray)
    for index,count in enumerate(frequencyArray):
        if(count==maxcount):
            result.append(numbertopattern(index,k))
    return result


def frequentwordswithmismatchessorted(text,k,d):
    frequentPatterns = []
    neighborhood = []
    index = []
    count = []
    for i in range(len(text)-k+1):
        neighborhood.append(Neighbors(text[i:i+k],d))
    neighborhood = [j for l1 in neighborhood for j in l1]
    #print(len(neighborhood))
    for i in range(len(neighborhood)):
        pattern = neighborhood[i]
        index.append(patterntonumber(pattern))
        count.append(1)
    sortedindex = index
    sortedindex.sort()
    for i in range(len(neighborhood)-1):
        if (sortedindex[i] == sortedindex[i+1]):
            count[i+1] = count[i]+1
    maxcount = max(count)
    for i in range(len(neighborhood)):
        if(count[i] == maxcount):
            pattern = numbertopattern(sortedindex[i], k)
            frequentPatterns.append(pattern)
    return(frequentPatterns)

l = frequentwordswithmismatchessorted('AAAATTCGCGATTATTCGCGCGCGCGCTCGCGAAAAAACGCGGAAAGAAAATTCGCTCGCTAAACGCGAAAAAACGCGATTATTAAAGAAACGCTATTGAAAAAAAAAAAAATTAAAAAACGCGCGCGAAAAAACGCGGAAAATTGAAACGCGGAAACGCGGAAACGCGCGCGCGCTGAAACGCGGAAACGCTCGCGAAAAAAATTCGCTCGCTATTATTGAAACGCTCGCGCGCTCGCGGAAACGCTGAAAATTCGCTGAAAGAAAAAAGAAAATTCGCGGAAAGAAACGCTAAACGCGCGCGATTATTATTAAACGCTCGCTCGCGAAAATTCGCGCGCTATTCGCTCGCTAAAATT',5,2)

def swap(letter):
    if letter == 'A':
        return 'T'
    elif letter == 'T':
        return 'A'
    elif letter == 'G':
        return 'C'
    else:
        return 'G'

def revcomplement(pattern):
    rev = pattern[::-1]
    comp=''
    for i in range(len(pattern)):
        comp = comp+swap(rev[i])
    return comp

def frequentwordswithmismatches_sorted_revcomp(text,k,d):
    frequentPatterns = []
    neighborhood = []
    index = []
    count = []
    for i in range(len(text)-k+1):
        neighborhood.append(Neighbors(text[i:i+k],d))
    neighborhood = [j for l1 in neighborhood for j in l1]
    #print(len(neighborhood))
    for i in range(len(neighborhood)):
        pattern = revcomplement(neighborhood[i])
        index.append(patterntonumber(pattern))
        count.append(1)
    for i in range(len(neighborhood)):
        pattern = neighborhood[i]
        index.append(patterntonumber(pattern))
        count.append(1)
    sortedindex = index
    sortedindex.sort()
    for i in range(len(neighborhood)-1):
        if (sortedindex[i] == sortedindex[i+1]):
            count[i+1] = count[i]+1
    maxcount = max(count)
    for i in range(len(neighborhood)):
        if(count[i] == maxcount):
            pattern = numbertopattern(sortedindex[i], k)
            frequentPatterns.append(pattern)
    return(frequentPatterns)

ans = frequentwordswithmismatches_sorted_revcomp('CACAAAGTATCAATCTGTGAAGCATGAAGAAGTGATCAAGATCCAAAGTATTATAAGTATATCTGAAGATCAAGCACAAAGAAGTGTGCACACATATTGCACATATCAATCATCTGAAGCATATATCATCTATATCTGATCAAGCATGATCATCCAATCTATAAGTGTATATCATCATCCATATAAGAAGTGAAGCACAATCTG',6,2)



with open('Salmonella_enterica.txt',"r") as file:
    text = file.read().split('\n')
    #print(text[0])
    genome=''
    for i in range(1,len(text)):
        genome+=text[i]
l,res = skew(genome)


def computingFrequencywithmissmatches(text,k,d):
    frequencyArray = []
    result = []
    neighborhood = []
    for i in range(4**k):
        frequencyArray.append(0)
    for i in range(len(text)-k+1):
        #pattern = text[i:i+k]
        neighborhood.append(Neighbors(text[i:i+k],d))
        #neighborhood = Neighbors(pattern,d)
    neighborhood = [j for l1 in neighborhood for j in l1]
    for str in neighborhood:
        j = patterntonumber(str)
        frequencyArray[j] +=1
    maxcount = max(frequencyArray)
    for index,count in enumerate(frequencyArray):
        if(count==maxcount):
            result.append(numbertopattern(index,k))
    return result


def computingFrequencywithmissmatchesrevcomp(text,k,d):
    frequencyArray = []
    result = []
    neighborhood = []
    for i in range(4**k):
        frequencyArray.append(0)
    for i in range(len(text)-k+1):
        #pattern = text[i:i+k]
        neighborhood.append(Neighbors(text[i:i+k],d))
        neighborhood.append(Neighbors(revcomplement(text[i:i+k]),d))
        #neighborhood = Neighbors(pattern,d)
    neighborhood = [j for l1 in neighborhood for j in l1]
    for str in neighborhood:
        j = patterntonumber(str)
        frequencyArray[j] +=1
    maxcount = max(frequencyArray)
    for index,count in enumerate(frequencyArray):
        if(count==maxcount):
            result.append(numbertopattern(index,k))
    return result
