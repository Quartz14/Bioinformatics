def patternCount(text, pattern):
    count = 0
    for i in range(len(text)-len(pattern)+1):
        if(text[i:i+len(pattern)] == pattern):
            count += 1
    return (count)

with open('dataset_2_7.txt',"r") as file:
    text,pattern,_ = file.read().split('\n')
    print(patternCount(text,pattern))

def FrequentWords(text, k):
    frequentPatterns = []
    count = []
    for i in range(len(text)-k+1):
        pattern = text[i:i+k]
        count.append(patternCount(text,pattern))
    maxcount = max(count)
    for i in range(len(text)-k+1):
        if (count[i] == maxcount):
            frequentPatterns.append(text[i:i+k])
    frequentPatterns = set(frequentPatterns)
    return(frequentPatterns)

with open("dataset_2_10.txt","r") as file:
    text,k,_ = file.read().split('\n')
    print(FrequentWords(text, int(k)))

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

with open("dataset_3_2.txt","r") as file:
    pattern,_ = file.read().split('\n')
    print(revcomplement(pattern)

def patmatch(pattern, text):
    match = []
    for i in range(len(text)-len(pattern)+1):
        if(text[i:i+len(pattern)] == pattern):
            match.append(i)
    return (match)
with open('dataset_3_5.txt',"r") as file:
    pattern,text,_ = file.read().split('\n')
    l = patmatch(pattern, text)
    print(*l) #Space delimited
with open('Vibrio_cholerae.txt',"r") as file:
    text,_ = file.read().split('\n')
    pattern = 'CTTGATCAT'
    l = patmatch(pattern, text)
    print(*l)



def FrequentWords_array(text, k):
    frequentPatterns = []
    count = []
    pattern_to_number = {}
    for i in range(len(text)-k+1):
        pattern = text[i:i+k]
        pattern_to_number[pattern] = pattern_to_number.get(pattern,0)+1
        #count.append(patternCount(text,pattern))
    maxcount = max(pattern_to_number.values())
    for key,val in pattern_to_number.items():
        if (val == maxcount):
            frequentPatterns.append(key)
    #frequentPatterns = set(frequentPatterns)
    return(frequentPatterns)

def symboltonumber(symbol):
    if (symbol=='A'):
        return(0)
    elif(symbol=='C'):
        return(1)
    elif(symbol=='G'):
        return(2)
    elif(symbol=='T'):
        return(3)

def numbertosymbol(number):
    if (number==0):
        return('A')
    elif(number==1):
        return('C')
    elif(number==2):
        return('G')
    elif(number==3):
        return('T')

def patterntonumber(pattern):
    if len(pattern) == 0:
        return (0)
    symbol = pattern[-1]
    prefix = pattern[:-1]
    return (4*patterntonumber(prefix) + symboltonumber(symbol))

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

def computingFrequency(text,k):
    frequencyArray = []
    for i in range(4**k):
        frequencyArray.append(0)
    for i in range(len(text)-k+1):
        pattern= text[i:i+k]
        j = patterntonumber(pattern)
        frequencyArray[j] = frequencyArray[j]+1
    return (frequencyArray)

def clumpFinding(genome,k,L,t):
    frequentPatterns = []
    clump = []
    for i in range(4**k):
        clump.append(0)
    for i in range(len(genome)-L+1):
        text = genome[i:i+L]
        frequencyArray = computingFrequency(text, k)
        for index in range(4**k):
            if(frequencyArray[index]>=t):
                #print(index,frequencyArray[index])
                clump[index] = 1
    for i in range(4**k):
        if clump[i] == 1:
            pattern = numbertopattern(i,k)
            #print(i,clump[i],pattern)
            frequentPatterns.append(pattern)
    return(frequentPatterns)

def BetterclumpFinding(genome,k,L,t):
    frequentPatterns = []
    clump = []
    for i in range(4**k):
        clump.append(0)
    text = genome[0:L]
    frequencyArray = computingFrequency(text, k)
    for i in range(4**k):
        if(frequencyArray[i]>=t):
                clump[i] = 1
    for i in range(1,len(genome)-L+1):
        firstpattern = genome[i-1:i-1+k]
        index = patterntonumber(firstpattern)
        frequencyArray[index] = frequencyArray[index]-1
        lastpattern = genome[i+L-k:i+L-k+k]
        index = patterntonumber(lastpattern)
        frequencyArray[index] = frequencyArray[index] + 1
        if frequencyArray[index] >= t:
            clump[index] = 1
    for i in range(4**k):
        if clump[i] == 1:
            pattern = numbertopattern(i,k)
            #print(i,clump[i],pattern)
            frequentPatterns.append(pattern)
    return(frequentPatterns)

with open('E_coli.txt',"r") as file:
    text = file.read()
    l = BetterclumpFinding(text,9,500,3)
    print(*l)
