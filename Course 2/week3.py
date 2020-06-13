
with open('RNA_codon_table_1.txt', 'r') as file:
    codon_protein_map = {}
    for line in file:
        lst = line.split()
        if(len(lst)>1):
            codon_protein_map[lst[0]] = lst[1]
        else:
            codon_protein_map[lst[0]] = ''


with open('integer_mass_table.txt', 'r') as file:
    protein_mass_map = {}
    for line in file:
        lst = line.split()
        protein_mass_map[lst[0]] = int(lst[1])


def protein_translate(genome):
    #print(genome)
    codons = [genome[i:i+3] for i in range(0,len(genome),3)]
    #print(codons)
    protein = ''
    for i in codons:
        protein += codon_protein_map[i]
    return protein

def mass(genome):
    m = 0
    for i in genome:
        m += protein_mass_map[i]
    return m

with open('dataset_96_4.txt','r') as file:
    text = file.read().strip()

protein_translate(text)

def translate(c):
    if(c == 'A'):
        return 'U'
    elif(c == 'T'):
        return 'A'
    elif(c == 'U'):
        return 'A'
    elif(c == 'G'):
        return 'C'
    elif (c == 'C'):
        return 'G'

def rev_translate_genome(genome):
    revgenome = genome[::-1]
    tgenome = ''
    for i in revgenome:
        id = translate(i)
        tgenome += id
    return(tgenome)


def peptide_encoding(text, peptide):
    np = len(peptide)
    text = text.replace('T','U')
    pep_dna = [text[i:i+(3*np)] for i in range(0,len(text))]
    filtered = [x for x in pep_dna if len(x) == 3*np]
    #print(filtered)
    l = []
    for p in filtered:
        cpeptide = protein_translate(p)
        rpeptide = protein_translate(rev_translate_genome(p))
        #print(p, cpeptide, rev_translate_genome(p),rpeptide)
        if((cpeptide == peptide) |(rpeptide == peptide)):
            p = p.replace('U','T')
            l.append(p)
    return(l)

with open('dataset_96_7.txt','r') as file:
with open('Bacillus_brevis.txt','r') as file:
    text = file.read().strip()
    text = text.replace('\n','')
    lst = text.split('\n')

cod = peptide_encoding(lst[0], lst[1])
cod = peptide_encoding(text,'VKLFPWFNQY')
print(len(cod))
for c in cod:
    print(c)

plst = peptide_encoding('ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA', 'MA')



protein_mass_map

def linearspectrum(peptide):
    prefixmass = [0]
    for i in range(0,len(peptide)):
        prefixmass.append(prefixmass[i] + protein_mass_map[peptide[i]])
    listspectrum = [0]#prefixmass
    for i in range(len(prefixmass)-1):
        for j in range(i+1, len(prefixmass)):
            listspectrum.append(prefixmass[j]-prefixmass[i])
    listspectrum.sort()
    return(listspectrum)

def cyclicspectrum(peptide):
    prefixmass = [0]
    for i in range(len(peptide)):
        prefixmass.append(prefixmass[i] + protein_mass_map[peptide[i]])
    peptidemass= prefixmass[-1]
    cyclicspectrum = [0]#prefixmass
    for i in range(len(prefixmass)-1):
        for j in range(i+1, len(prefixmass)):
            cyclicspectrum.append(prefixmass[j]-prefixmass[i])
            if((i>0) & (j<len(prefixmass)-1)):
                cyclicspectrum.append(peptidemass - (prefixmass[j] - prefixmass[i]))
    cyclicspectrum.sort()
    return(cyclicspectrum)

def consistent(peptide, spectrum):
    l = linearspectrum(peptide)
    dict = {}
    sdict= {}
    for i in l:
        dict[i] = dict.get(i,0)+1
    #print(dict.items())
    for i in spectrum:
    #    print(i)
        sdict[i] = sdict.get(i,0)+1
    #print(sdict.keys())
    for k,v in dict.items():
    #    print(k,v)
        if (k not in list(sdict.keys())):
            #print('not there')
            #print(k)
            return False
        elif(v > sdict[k]):
            #print('not enough')
            #print(v,sdict[k])
            return False
    return True

def expand(peptides):
    expanded_pep = []
    if(len(peptides)>0):
        for i in peptides:
            for j in ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'N', 'D', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']:
                expanded_pep.append(i+j)
    else:
        return(['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'N', 'D', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W'])
    return expanded_pep





d = [0 ,97, 99 ,113 ,114 ,128 ,128, 147, 147, 163, 186, 227, 241, 242, 244 ,260, 261, 262, 283, 291 ,333, 340, 357, 388 ,389 ,390 ,390 ,405 ,430, 430, 447, 485, 487, 503, 504, 518, 543, 544 ,552 ,575 ,577 ,584, 631, 632, 650, 651 ,671 ,672 ,690 ,691 ,738 ,745 ,747 ,770 ,778 ,779 ,804 ,818 ,819 ,835, 837, 875 ,892, 892, 917, 932, 932, 933 ,934 ,965 ,982 ,989 ,1031 ,1039 ,1060 ,1061 ,1062 ,1078 ,1080 ,1081 ,1095 ,1136 ,1159 ,1175, 1175 ,1194 ,1194, 1208 ,1209 ,1223 ,1225 ,1322]


s = '0 71 101 113 131 184 202 214 232 285 303 315 345 416'
sd = s.split()
sd = [int(i) for i in sd]

def cyclopeptidesequencing(spectrum):
    candidatepeptides = []
    flag = 0
    eflag = 1
    finalpeptides = []
    while ((len(candidatepeptides)>0) | (flag==0)):
        temp = []
        #if(eflag):
        candidatepeptides = expand(candidatepeptides)
        flag = 1
        for peptide in candidatepeptides:
            if(mass(peptide) == max(spectrum)):
                if((cyclicspectrum(peptide) == spectrum) & (peptide not in finalpeptides)):
                    finalpeptides.append(peptide)
                    #print(peptide)
                #candidatepeptides.remove(peptide)
                #eflag = 0
            #elif not (consistent(peptide,spectrum)):
                #candidatepeptides.remove(peptide)
            if(consistent(peptide,spectrum)):
                temp.append(peptide)
        candidatepeptides = list(temp)
        #print(candidatepeptides)
    return(finalpeptides)
fp = cyclopeptidesequencing(text)
ans = []
for i in fp:
    s = '-'.join([str(protein_mass_map[l]) for l in i])
    ans.append(s)


print(*ans)

with open('leaderboard_spectrum.txt','r') as file:
with open('dataset_100_6.txt','r') as file:
    text = file.read().strip()
    text = text.split(' ')
    text = [int(i) for i in text]



Alphabet = {57: 'G', 71: 'A', 87: 'S', 97: 'P',
           99: 'V', 101: 'T', 103: 'C', 113:'I/L',
           114: 'N', 115: 'D', 128: 'K/Q', 129: 'E',
           131: 'M', 137: 'H', 147: 'F', 156: 'R', 163: 'Y', 186: 'W'}

def CountingMass(Mass, masslist):
    if Mass == 0: return 1, masslist
    if Mass < 57: return 0, masslist
    if Mass in masslist: return masslist[Mass], masslist
    n = 0
    for i in Alphabet:
        k, masslist = CountingMass(Mass - i, masslist)
        n += k
    masslist[Mass] = n
    return n, masslist

print(CountingMass(1024, {})[0])
