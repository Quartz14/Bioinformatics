def cyscoring(peptide, spectrum):
    #l = linearspectrum(peptide)
    l = cyclicspectrum(peptide)
    score = 0
    dict = {}
    sdict= {}
    for i in l:
        dict[i] = dict.get(i,0)+1
    for i in spectrum:
        sdict[i] = sdict.get(i,0)+1
    for k,v in dict.items():
        if (k in list(sdict.keys())):
            score += min([v,sdict[k]])
    return score

def linearscoring(peptide, spectrum):
    l = linearspectrum(peptide)
    score = 0
    dict = {}
    sdict= {}
    for i in l:
        dict[i] = dict.get(i,0)+1
    for i in spectrum:
        sdict[i] = sdict.get(i,0)+1
    for k,v in dict.items():
        if (k in list(sdict.keys())):
            score += min([v,sdict[k]])
    return score

spectrum = [0 ,99 ,113 ,114 ,128 ,227 ,257 ,299 ,355 ,356 ,370 ,371, 484]
pep = 'NQEL'
cyscoring(pep,spectrum)

with open('dataset_4913_1.txt','r') as file:
    text = file.read().strip()
    text = text.split('\n')
    pep = text[0]
    spec = [int(i) for i in text[1].split(' ')]

with open('dataset_102_3.txt','r') as file:
    text = file.read().strip()
    text = text.split('\n')
    pep = text[0]
    spec = [int(i) for i in text[1].split(' ')]

def leaderboardcyclopeptidesequencing(spectrum, n):
    leaderboard = []
    leaderpeptide = ''
    lp = {0:['']}
    flag = 1
    while((len(leaderboard)>0) | flag ==1):
        flag=0
        leaderboard = expand(leaderboard)
        dummylb = list(leaderboard)
        for peptide in leaderboard:
            #print(peptide)
            if(mass(peptide) == spectrum[-1]):
                #if(linearscoring(peptide, spectrum) > linearscoring(leaderpeptide, spectrum)):
                if(cyscoring(peptide, spectrum) >= cyscoring(leaderpeptide, spectrum)):
                    leaderpeptide = peptide
                    #print(peptide,cyscoring(peptide, spectrum))
                    #print(lp)
                    if(cyscoring(peptide, spectrum) not in lp.keys()):
                        lp[cyscoring(peptide, spectrum)] = []
                    lp[cyscoring(peptide,spectrum)].append(leaderpeptide)
                    #print(peptide,cyscoring(peptide, spectrum))
            elif (mass(peptide)>spectrum[-1]):
                dummylb.remove(peptide)
        #print(leaderboard)
        #print(dummylb)
        leaderboard = trim(dummylb,spectrum,n)
        #print(leaderboard)
    return lp[max(lp)]#leaderpeptide

with open('dataset_102_10.txt','r') as file:
    text = file.read().strip()
    text = text.split('\n')
    n = int(text[0])
    spec = [int(i) for i in text[1].split(' ')]

spec = [0 ,71 ,113, 129 ,147 ,200, 218 ,260 ,313 ,331 ,347 ,389 ,460]
n  =10
with open('dataset_102_8.txt','r') as file:
    text = file.read().strip()
    text = text.split('\n')
    n = int(text[0])
    spec = [int(i) for i in text[1].split(' ')]

lp = leaderboardcyclopeptidesequencing(spec, n)
ans = []
for i in lp:
    s = '-'.join([str(protein_mass_map[l]) for l in i])
    ans.append(s)


print(*ans)

with open('Tyrocidine_B1_Spectrum_10.txt','r') as file:
with open('Tyrocidine_B1_Spectrum_25.txt','r') as file:
    text = file.read().strip()
    spec = [int(i) for i in text.split(' ')]



def trim(leaderboard, spectrum, n):
    linearscores = []
    for i in leaderboard:
        linearscores.append(0)
    for j in range(len(leaderboard)):
        peptide = leaderboard[j]
        linearscores[j] = linearscoring(peptide, spectrum)
    sleaderboard = [x for _,x in sorted(zip(linearscores,leaderboard),reverse=True)]
    #print(sleaderboard)
    linearscores.sort(reverse=True)
    #print(linearscores)
    for j in range(n+1,len(leaderboard)):
        if(linearscores[j]< linearscores[n]):
            return(sleaderboard[:j-1])
    return sleaderboard

lb = ['LAST' ,'ALST' ,'TLLT' ,'TQAS']
spec = [0 ,71 ,87 ,101 ,113 ,158 ,184 ,188 ,259 ,271 ,372]
n = 2
nlb = trim(lb,spec,n)
nlb

with open('dataset_4913_3.txt','r') as file:
    text = file.read().strip()
    text = text.split('\n')
    leaderboard = [i for i in text[0].split(' ')]
    spec = [int(i) for i in text[1].split(' ')]
    n = int(text[2])

lb = trim(leaderboard, spec, n)
print(*lb)


def convolution(spectrum):
    l = []
    spectrum.sort()
    for i in range(len(spectrum)-1):
        for j in range(i+1,len(spectrum)):
            #print(spectrum[j],spectrum[i],spectrum[j]-spectrum[i])
            l.append(spectrum[j]-spectrum[i])
    #l.sort()
    lcopy = list(l)
    for i in l:
        if(i == 0):
            lcopy.remove(i)
    return(lcopy)


lconv = convolution([0,137,186,323])
print(*lconv)

with open('dataset_104_4.txt','r') as file:
    text = file.read().strip()
    spec = [int(i) for i in text.split(' ')]

lconv = convolution(spec)
print(*lconv)

def convolutioncyclopeptidesequencing(m,n,spectrum):
    lst = convolution(spectrum)
    lmax = max(spectrum)
    #print(lmax)
    #print('--------------------------------')
    newlst = []
    for i in lst:
        if((i>=57) & (i<200)):
            newlst.append(i)
    #print(newlst)
    #print('--------------------------------')
    newlst1 = collections.Counter(newlst).most_common(m)
    #print(newlst1)
    #print('--------------------------------')
    frequency = newlst1[-1][1]
    #print(frequency)
    #print('--------------------------------')
    newlst2 = dict(collections.Counter(newlst))
    #print(newlst2)
    #print('--------------------------------')
    extendedalphabetmass = [k for k,v in newlst2.items() if v >= frequency]
    extendedalphabetmass.sort()
    #print(extendedalphabetmass)
    #print('--------------------------------')
    #print(len(extendedalphabetmass))
    #print('--------------------------------')
    myextended_list = {}
    for i in extendedalphabetmass:
        myextended_list[chr(i)] = i
    print(myextended_list)
    print('--------------------------------')
    ans = leaderboardcyclopeptidesequencing_extended(extendedalphabetmass , n,myextended_list, lmax)
    return (ans, myextended_list)
ll = [57 ,57 ,71 ,99 ,129 ,137 ,170 ,186 ,194 ,208 ,228 ,265 ,285 ,299 ,307 ,323 ,356, 364 ,394 ,422 ,493]
lspec, extended_list = convolutioncyclopeptidesequencing(20,60,ll)

with open('dataset_104_7.txt','r') as file:
    text = file.read().strip()
    text = text.split('\n')
    m = int(text[0])
    n = int(text[1])
    spec = [int(i) for i in text[2].split(' ')]

lspec, extended_list = convolutioncyclopeptidesequencing(m,n,spec)
ans = []
for i in lspec:
    s = '-'.join([str(extended_list[l]) for l in i])
    ans.append(s)


print(*ans)


aminoacid = []
aminoacidmass = []
for i in range(57,201):
    aminoacidmass.append(i)
    aminoacid += chr(i)


def expand_extended(peptides,extended_list=extended_list):
    expanded_pep = []
    if(len(peptides)>0):
        for i in peptides:
            for j in extended_list.keys():
                expanded_pep.append(i+j)
    else:
        return(extended_list.keys())
    return expanded_pep

extended_list = {}
for i in aminoacidmass:
    extended_list[chr(i)] = i

def mass_extended(genome,extended_list=extended_list):
    m = 0
    for i in genome:
        m += extended_list[i]
    return m
def cyclicspectrum_extended(peptide,extended_list=extended_list):
    prefixmass = [0]
    for i in range(len(peptide)):
        prefixmass.append(prefixmass[i] + extended_list[peptide[i]])
    peptidemass= prefixmass[-1]
    cyclicspectrum = [0]#prefixmass
    for i in range(len(prefixmass)-1):
        for j in range(i+1, len(prefixmass)):
            cyclicspectrum.append(prefixmass[j]-prefixmass[i])
            if((i>0) & (j<len(prefixmass)-1)):
                cyclicspectrum.append(peptidemass - (prefixmass[j] - prefixmass[i]))
    cyclicspectrum.sort()
    return(cyclicspectrum)
def cyscoring_extended(peptide, spectrum,extended_list=extended_list):
    #l = linearspectrum(peptide)
    l = cyclicspectrum_extended(peptide,extended_list=extended_list)
    score = 0
    dict = {}
    sdict= {}
    for i in l:
        dict[i] = dict.get(i,0)+1
    for i in spectrum:
        sdict[i] = sdict.get(i,0)+1
    for k,v in dict.items():
        if (k in list(sdict.keys())):
            score += min([v,sdict[k]])
    return score

def linearspectrum_extended(peptide,extended_list=extended_list):
    prefixmass = [0]
    for i in range(0,len(peptide)):
        prefixmass.append(prefixmass[i] + extended_list[peptide[i]])
    listspectrum = [0]#prefixmass
    for i in range(len(prefixmass)-1):
        for j in range(i+1, len(prefixmass)):
            listspectrum.append(prefixmass[j]-prefixmass[i])
    listspectrum.sort()
    return(listspectrum)

def linearscoring_extended(peptide, spectrum,extended_list=extended_list):
    l = linearspectrum_extended(peptide,extended_list=extended_list)
    score = 0
    dict = {}
    sdict= {}
    for i in l:
        dict[i] = dict.get(i,0)+1
    for i in spectrum:
        sdict[i] = sdict.get(i,0)+1
    for k,v in dict.items():
        if (k in list(sdict.keys())):
            score += min([v,sdict[k]])
    return score

def trim_extended(leaderboard, spectrum, n,extended_list=extended_list):
    linearscores = []
    for i in leaderboard:
        linearscores.append(0)
    for j in range(len(leaderboard)):
        peptide = leaderboard[j]
        linearscores[j] =linearscoring_extended(peptide, spectrum,extended_list=extended_list)
        #print('trimmed')
        #print(peptide, linearscores[j])
    sleaderboard = [x for _,x in sorted(zip(linearscores,leaderboard),reverse=True)]
    #print(sleaderboard)
    linearscores.sort(reverse=True)
    #print(linearscores)
    for j in range(n+1,len(leaderboard)):
        if(linearscores[j]< linearscores[n]):
            return(sleaderboard[:j-1])
    return sleaderboard

def leaderboardcyclopeptidesequencing_extended(spectrum, n, extended_list, lmax):
    #print(f'-------------------{lmax}----------------------------')
    leaderboard = []
    leaderpeptide = ''
    lp = {0:['']}
    flag = 1
    while((len(leaderboard)>0) | flag ==1):
        flag=0
        leaderboard = expand_extended(leaderboard,extended_list=extended_list)
        #print(leaderboard)
        dummylb = list(leaderboard)
        for peptide in leaderboard:
            #print(peptide)
            if(mass_extended(peptide,extended_list=extended_list) == lmax):
                #if(linearscoring(peptide, spectrum) > linearscoring(leaderpeptide, spectrum)):
                #print(peptide,mass_extended(peptide,extended_list=extended_list))
                #print("************************************")
                if(cyscoring_extended(peptide, spectrum,extended_list=extended_list) >= cyscoring_extended(leaderpeptide, spectrum,extended_list=extended_list)):
                    leaderpeptide = peptide
                    #print(peptide,cyscoring(peptide, spectrum))
                    #print(leaderboard)
                    #print(leaderpeptide,cyscoring_extended(peptide, spectrum,extended_list=extended_list))
                    if(cyscoring_extended(peptide, spectrum,extended_list=extended_list) not in lp.keys()):
                        lp[cyscoring_extended(peptide, spectrum,extended_list=extended_list)] = []
                    lp[cyscoring_extended(peptide,spectrum,extended_list=extended_list)].append(leaderpeptide)
                    #print(peptide,cyscoring(peptide, spectrum))
            elif (mass_extended(peptide,extended_list=extended_list)>lmax):
                dummylb.remove(peptide)
        #print(leaderboard)
        #print(dummylb)
        leaderboard = trim_extended(dummylb,spectrum,n,extended_list=extended_list)
        #print(leaderboard)
    return lp#(lp[max(lp)])#leaderpeptide


with open('Tyrocidine_B1_Spectrum_10.txt','r') as file:
    text = file.read().strip()
    spec = [int(i) for i in text.split(' ')]

lpp = leaderboardcyclopeptidesequencing_extended(spec, n=1000)
ans = []
for i in lpp:
    s = '-'.join([str(extended_list[l]) for l in i])
    ans.append(s)
