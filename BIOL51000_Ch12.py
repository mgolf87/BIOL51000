###############################################################################
#####            CHAPTER 12: PAIRWISE SEQUENCE ALIGNMENTS                 #####
###############################################################################

#   Why Align?
#       1) evolutionary relationship - common ancestor
#       2) conserved biological function - degree of similarity
#       3) starting point to guide experiments

# Sequence classification = sequence annotation
# Natural divergence of sequences through DNA replication and crossing over
# Consider mutations as the tolerated changes - non-tolerated changes likely
# result in the termination or death of cell/indiv - look for non-variable genes
# in diseases instead of 'strong effect mutations'

# Protein Folding:
# active site = conserved, inner hydrophobic core = conserved, residues on surface
# and flexible regions = less conserved, as long as change doesnt fuck 3D shape

# Human chromosome 2 result of merging of two smaller chromosomes
# Sequences linked by ancestry = family, generally arise from sections of 
# chromosomes duplicated
# Genes which don't follow main evolutionary trend = rogue genes, result of
# genetic material trasmitted other method than reproduction - virus or symbiosis


###############################################################################
#####                        SEQUENCE IDENTITY                            #####
###############################################################################

def calcSeqIdentity (seqA, seqB):
    """ numPlaces: length of shortest sequence - guarantee we don't overshoot
        the smallest of the pair"""
        
    numPlaces = min(len(seqA), len(seqB))
    score = 0.0
        
    for i in range(numPlaces):
        if seqA[i] == seqB[i]:
            score += 1.0
        
    return 100.0* score/numPlaces

# Test:
seq1 = 'ALIGNMENTS'
seq2 = 'ALIGDVENTS'
seq3 = 'ALIGDPVENTS'
seq4 = 'ALIGN-MENTS'

print(calcSeqIdentity(seq1, seq2)) # 80.0
print(calcSeqIdentity(seq1, seq3)) # 40.0 
print(calcSeqIdentity(seq4, seq3)) # 72.7272
# illustrates how inserting a gap in the right place is crucial for alignment


###############################################################################
#####                        SUBSTITUTABILITY                             #####
###############################################################################

# residue substitution scores stored in 2D array = substitution/similarity matrix
# matrices = dictionaries of dictionaries -> residue letters directly as keys
# Substition matrix dictionary first key (residueletter) = sub-dictionary inside
# main dictionary, second key = final value within the sub-dictionary

DNA_1 = {'G': { 'G':1, 'C':0, 'A':0, 'T':0 },
         'C': { 'G':0, 'C':1, 'A':0, 'T':0 },
         'A': { 'G':0, 'C':0, 'A':1, 'T':0 },
         'T': { 'G':0, 'C':0, 'A':0, 'T':1 }}
DNA_1['G']['G'] # =1 as identical
DNA_1['G']['A'] # =0 as not

# Complementarity score, 1 for A:T or G:C and -1 for mismatches
REV_COMP = {'G': { 'G':-1, 'C': 1, 'A':-1, 'T':-1 },
            'C': { 'G': 1, 'C':-1, 'A':-1, 'T':-1 },
            'A': { 'G':-1, 'C':-1, 'A':-1, 'T': 1 },
            'T': { 'G':-1, 'C':-1, 'A': 1, 'T':-1 }}

# Substituion scores where N is unknown residue
DNA_2 = {'G': { 'G': 8, 'C':-5, 'A':-5, 'T':-5, 'N':0 },
         'C': { 'G':-5, 'C': 8, 'A':-5, 'T':-5, 'N':0 },
         'A': { 'G':-5, 'C':-5, 'A': 8, 'T':-5, 'N':0 },
         'T': { 'G':-5, 'C':-5, 'A':-5, 'T': 8, 'N':0 },
         'N': { 'G': 0, 'C': 0, 'A': 0, 'T': 0, 'N':0 }}  

# Matrix is symmetric BLOSUM62['A']['R'] = BLOSUM62['R']['A']
BLOSUM62 = {'A':{'A': 4,'R':-1,'N':-2,'D':-2,'C': 0,'Q':-1,'E':-1,'G': 0,'H':-2,'I':-1,
                 'L':-1,'K':-1,'M':-1,'F':-2,'P':-1,'S': 1,'T': 0,'W':-3,'Y':-2,'V': 0,'X':0},
            'R':{'A':-1,'R': 5,'N': 0,'D':-2,'C':-3,'Q': 1,'E': 0,'G':-2,'H': 0,'I':-3,
                 'L':-2,'K': 2,'M':-1,'F':-3,'P':-2,'S':-1,'T':-1,'W':-3,'Y':-2,'V':-3,'X':0},
            'N':{'A':-2,'R': 0,'N': 6,'D': 1,'C':-3,'Q': 0,'E': 0,'G': 0,'H': 1,'I':-3,
                 'L':-3,'K': 0,'M':-2,'F':-3,'P':-2,'S': 1,'T': 0,'W':-4,'Y':-2,'V':-3,'X':0},
            'D':{'A':-2,'R':-2,'N': 1,'D': 6,'C':-3,'Q': 0,'E': 2,'G':-1,'H':-1,'I':-3,
                 'L':-4,'K':-1,'M':-3,'F':-3,'P':-1,'S': 0,'T':-1,'W':-4,'Y':-3,'V':-3,'X':0},
            'C':{'A': 0,'R':-3,'N':-3,'D':-3,'C': 9,'Q':-3,'E':-4,'G':-3,'H':-3,'I':-1,
                 'L':-1,'K':-3,'M':-1,'F':-2,'P':-3,'S':-1,'T':-1,'W':-2,'Y':-2,'V':-1,'X':0},
            'Q':{'A':-1,'R': 1,'N': 0,'D': 0,'C':-3,'Q': 5,'E': 2,'G':-2,'H': 0,'I':-3,
                 'L':-2,'K': 1,'M': 0,'F':-3,'P':-1,'S': 0,'T':-1,'W':-2,'Y':-1,'V':-2,'X':0},
            'E':{'A':-1,'R': 0,'N': 0,'D': 2,'C':-4,'Q': 2,'E': 5,'G':-2,'H': 0,'I':-3,
                 'L':-3,'K': 1,'M':-2,'F':-3,'P':-1,'S': 0,'T':-1,'W':-3,'Y':-2,'V':-2,'X':0},
            'G':{'A': 0,'R':-2,'N': 0,'D':-1,'C':-3,'Q':-2,'E':-2,'G': 6,'H':-2,'I':-4,
                 'L':-4,'K':-2,'M':-3,'F':-3,'P':-2,'S': 0,'T':-2,'W':-2,'Y':-3,'V':-3,'X':0},
            'H':{'A':-2,'R': 0,'N': 1,'D':-1,'C':-3,'Q': 0,'E': 0,'G':-2,'H': 8,'I':-3,
                 'L':-3,'K':-1,'M':-2,'F':-1,'P':-2,'S':-1,'T':-2,'W':-2,'Y': 2,'V':-3,'X':0},
            'I':{'A':-1,'R':-3,'N':-3,'D':-3,'C':-1,'Q':-3,'E':-3,'G':-4,'H':-3,'I': 4,
                 'L': 2,'K':-3,'M': 1,'F': 0,'P':-3,'S':-2,'T':-1,'W':-3,'Y':-1,'V': 3,'X':0},
            'L':{'A':-1,'R':-2,'N':-3,'D':-4,'C':-1,'Q':-2,'E':-3,'G':-4,'H':-3,'I': 2,
                 'L': 4,'K':-2,'M': 2,'F': 0,'P':-3,'S':-2,'T':-1,'W':-2,'Y':-1,'V': 1,'X':0},
            'K':{'A':-1,'R': 2,'N': 0,'D':-1,'C':-3,'Q': 1,'E': 1,'G':-2,'H':-1,'I':-3,
                 'L':-2,'K': 5,'M':-1,'F':-3,'P':-1,'S': 0,'T':-1,'W':-3,'Y':-2,'V':-2,'X':0},
            'M':{'A':-1,'R':-1,'N':-2,'D':-3,'C':-1,'Q': 0,'E':-2,'G':-3,'H':-2,'I': 1,
                 'L': 2,'K':-1,'M': 5,'F': 0,'P':-2,'S':-1,'T':-1,'W':-1,'Y':-1,'V': 1,'X':0},
            'F':{'A':-2,'R':-3,'N':-3,'D':-3,'C':-2,'Q':-3,'E':-3,'G':-3,'H':-1,'I': 0,
                 'L': 0,'K':-3,'M': 0,'F': 6,'P':-4,'S':-2,'T':-2,'W': 1,'Y': 3,'V':-1,'X':0},
            'P':{'A':-1,'R':-2,'N':-2,'D':-1,'C':-3,'Q':-1,'E':-1,'G':-2,'H':-2,'I':-3,
                 'L':-3,'K':-1,'M':-2,'F':-4,'P': 7,'S':-1,'T':-1,'W':-4,'Y':-3,'V':-2,'X':0},
            'S':{'A': 1,'R':-1,'N': 1,'D': 0,'C':-1,'Q': 0,'E': 0,'G': 0,'H':-1,'I':-2,
                 'L':-2,'K': 0,'M':-1,'F':-2,'P':-1,'S': 4,'T': 1,'W':-3,'Y':-2,'V':-2,'X':0},
            'T':{'A': 0,'R':-1,'N': 0,'D':-1,'C':-1,'Q':-1,'E':-1,'G':-2,'H':-2,'I':-1,
                 'L':-1,'K':-1,'M':-1,'F':-2,'P':-1,'S': 1,'T': 5,'W':-2,'Y':-2,'V': 0,'X':0},
            'W':{'A':-3,'R':-3,'N':-4,'D':-4,'C':-2,'Q':-2,'E':-3,'G':-2,'H':-2,'I':-3,
                 'L':-2,'K':-3,'M':-1,'F': 1,'P':-4,'S':-3,'T':-2,'W':11,'Y': 2,'V':-3,'X':0},
            'Y':{'A':-2,'R':-2,'N':-2,'D':-3,'C':-2,'Q':-1,'E':-2,'G':-3,'H': 2,'I':-1,
                 'L':-1,'K':-2,'M':-1,'F': 3,'P':-3,'S':-2,'T':-2,'W': 2,'Y': 7,'V':-1,'X':0},
            'V':{'A': 0,'R':-3,'N':-3,'D':-3,'C':-1,'Q':-2,'E':-2,'G':-3,'H':-3,'I': 3,
                 'L': 1,'K':-2,'M': 1,'F':-1,'P':-2,'S':-2,'T': 0,'W':-3,'Y':-1,'V': 4,'X':0},
            'X':{'A': 0,'R': 0,'N': 0,'D': 0,'C': 0,'Q': 0,'E': 0,'G': 0,'H': 0,'I': 0,
                 'L': 0,'K': 0,'M': 0,'F': 0,'P': 0,'S': 0,'T': 0,'W': 0,'Y': 0,'V': 0,'X':0}}
# Diagonal is not uniform, ['A']['A'] = 4, ['N']['N'] = 6
# 'A' is less well conserved = more swappable than 'N'

###############################################################################
#####                        SIMILARITY                                   #####
###############################################################################

def calcSeqSimilarity(seqA, seqB, simMatrix):
    
    numPlaces = min(len(seqA), len(seqB))
    totalScore = 0.0
    
    for i in range(numPlaces):
        residueA = seqA[i]
        residueB = seqB[i]
        
        totalScore += simMatrix[residueA][residueB]
    return totalScore
# looks up a similarity score in substitution matrix rather than checking
# equivalence, cannot deal with '-' gaps though
# Test:
# DNA
print(calcSeqSimilarity('AGCATCGCTCT', 'AGCATCGTTTT', DNA_2)) # 62.0
# PROTEIN
print(calcSeqSimilarity('ALIGNMENT', 'AYIPNVENT', BLOSUM62)) # 28.0

def pairAlignScore (alignA, alignB, simMatrix, insert=8, extend=4):
    totalScore = 0.0
    n = min(len(alignA), len(alignB))
    
    for i in range(n):
        residueA = alignA[i]
        residueB = alignB[i]
        
        if '-' not in (residueA, residueB): # not a gap
            simScore = simMatrix[residueA][residueB]
        elif (i > 0) and ('-' in (alignA[i-1], alignB[i-1])): # if gap not at start and previous position a gap we substract extend penalty from total score
            simScore = -extend
        else: # new gap if either condition not met above substracting the insert penalty
            simScore = -insert
        
        totalScore += simScore
    return totalScore

# Test
print(pairAlignScore('ALIGDPPVENTS', 'ALIGN--MENTS', BLOSUM62)) # 28.0
print(pairAlignScore('ALIGDPPVENTS', '--ALIGNMENTS', BLOSUM62)) # -3.0

###############################################################################
#####                        PAIRWISE ALIGNMENT                           #####
###############################################################################

# Optimising - which alginment out of all possible combinations is the best (highest score)
# Same alignments represented as a comparison matrix - each element represents
# the pairing of a residue of one sequence with another
# aligned (paired) residues for given row and column with 'x'
# alternative alignments are different routes within the matrix, gaps are jumps
# down rows or across columns

### DYNAMIC PROGRAMMING - breaking a big problem into sub-problems which can be solved
# once and are often repeated within the larger problem
# sub-problems only follow best paths not all paths (kind of quantum-like)
# goal: route of maximum alignment score (not min time)
# three routes into each point: (1) gap in one sequence, (2) gap in the other 
# sequence, (3) align two residues
# repeatedly extending alignment and pruning sub-optimal alternatives
# We work backwards from the best alignment because as the sequences are compared
# we dont know which sub-alignments are discarded or which win out but we know
# what decisions are made as they are recorded

def sequenceAlign(seqA, seqB, simMatrix=DNA_2, insert=8, extend=4):
    numI = len(seqA) + 1 # represent the size of the comparison grid
    numJ = len(seqB) + 1 # one larger as require extra row column for starting values
    
    scoreMatrix = [[0] * numJ for x in range(numI)] # lists of list but could use NumPy
    routeMatrix = [[0] * numJ for x in range(numI)] # 0 = pairing, 1 = gap in seqB, 2 = gap in seqA

    for i in range(1, numI): # adjust top and left edges except first element
        routeMatrix[i][0] = 1 # route codes = gaps, good for indented sequences
    for j in range(1, numJ):
        routeMatrix[0][j] = 2
        
    for i in range(1, numI): # fill remainder of matrix
        for j in range(1, numJ):
            penalty1 = insert
            penalty2 = insert
            
            if routeMatrix[i-1][j] == 1: # detects whether previous position has a gap
                penalty1 = extend
            elif routeMatrix[i][j-1] == 2:
                penalty2 = extend
            similarity = simMatrix[seqA[i-1]][seqB[j-1]] # position for residues within two sequences - 1 larger due to extra initialisation values at start
            
            paths = [scoreMatrix[i-1][j-1] + similarity, # Route 0 - similarity score of two residues - going diagonally in comparison matrix i-1, j-1 to i,j
                     scoreMatrix[i-1][j] - penalty1, # Route 1 - gap in seqB - goind down a row i-1, j to i,j
                     scoreMatrix[i][j-1] - penalty2] # Route 2 - gap in seqA - across column i, j-1 to i,j
            best = max(paths) # max value in paths list 
            route = paths.index(best) # route code is index of this score, paths.index is the reason for using numeric codes
            
            scoreMatrix[i][j] = best
            routeMatrix[i][j] = route # updated as it iterates through
    
    alignA = [] # follow winning routes backwards - each point adding appropriate
    alignB = [] # gap or paired residues - each input sequecne filled with residue codes and dashes in reverse order
    i = numI-1 # row and column for the end of the alignment
    j = numJ-1
    score = scoreMatrix[i][j]
    
    while i > 0 or j > 0: # fills in alignment strings by going back along rows/columns
        route = routeMatrix[i][j] # only stops when both are 0 i.e. gaps
        if route == 0: # Diagonal
            alignA.append(seqA[i-j])
            alignB.append(seqB[j-1])
            i -= 1
            j -= 1
        elif route == 1: # Gap in seqB
            alignA.append(seqA[i-1])
            alignB.append('-')
            i -= 1
        elif route == 2: # Gap in seqA
            alignA.append('-')
            alignB.append(seqB[j-1])
            j -= 1
    
    alignA.reverse() # reverse alignment lists - as we were working backwards
    alignB.reverse()
    alignA = ''.join(alignA)
    alignB = ''.join(alignB)
    
    return score, alignA, alignB
# TEST:
seqA = 'WFSEPEIST'
seqB = 'FSRPAVVIST'

score, alignA, alignB = sequenceAlign(seqA, seqB, BLOSUM62)
print(score) # 17
print(alignA) # WFSEPE--IST
print(alignB) # -FSRPAVVIST

###############################################################################
#####                              BLAST                                  #####
###############################################################################

from re import sub # formats sequnces for input files
from subprocess import call # runs programs external to Python
from xml.etree import ElementTree # let's us read XML formatted BLAST output
# used in the second example listed but should be contained within the code files
def makeBlastDatabase(fastaFile, databaseName, 
                      formatdbExe=None, isProtein=True):
    """ specially formatted version of sequence data to be queried, which allows
        BLAST to efficiently find matches. 
        FASTA-formatted file as input
        formatdbExe = database
        location of external 'makeblastdb' program to run and 
        isProtein --> False = DNA/RNA """
        # formatdbExe value not passed -> assume database creation available to local system as command makeblastdb
    if not formatdbExe:
        formatdbExe = 'makeblastdb'
    if isProtein: # defines variable passed and makeblastdb accepts below
        molType = 'prot'
    else:
        molType = 'nuc1'
        
    cmdArgs = [formatdbExe,
               '-dbtype', molType,
               '-in', fastaFile,
               '-out', databaseName] # list passed to call function to run
    print('Making BLAST database %s...' % databaseName)
    
    try:
        call(cmdArgs)
    except Exception as err: # exception triggered - print command tried and error
        print('BLAST database creation failed')
        print('Command used: "%s"'%''.join(cmdArgs))
        print(err)
        return
    print('...done')
    # defines function, can run makeblastdb program on FASTA format file - must be done once per database
    
    fileName = 'EcoliGenome.fasta'
    makeBlastDatabase(fileName, 'ECOLI_PROT')
    #makeBlastDatabase(fileName, 'ECOLI_PROT', '/usr/bin/makeblastdb') #if system doesnt know location of executable function
    
        












































