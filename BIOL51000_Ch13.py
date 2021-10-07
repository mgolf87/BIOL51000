###############################################################################
#####               CHAPTER 13: MULTIPLE-SEQUENCE ALIGNMENTS              #####
###############################################################################

# Closely related sequences: differences are significant, similarities are expected
# Distantly related sequences: differences are expected, similarities are significant
# Pairwise alignment = 2D - Dynamic Programming - every other sequence added
# adds another dimension of complexity

# PROGRESSIVE PAIRING: join sequences into alignment in order of similarity
# done via - construction of family tree for input sequences leading to creation
# of overall multiple alignment via progressively pairing alignments of smaller
# numbers of sequences in the same order as branches of family tree join#
# weighted dependent on the lengths of tree branchs - PHYLOGENETIC TREES
# Search speed-up: when doing exhaustive search restrict for optimal solutions
# so only routes lie near initial quick solution


###############################################################################
#####                             CONSENSUS                               #####
###############################################################################
# from multiple alignments we will produce a consensus sequence
# ASSUMPTIONS: input alignment list of one-letter sequences with gaps, each list
# assumed to be of the same length, threshold value input = position undefined
# if not met code 'X' for residue - default 0.25

def consensus(alignment, threshold=0.25):
    
    n = len(alignment[0]) # first sequence is measured
    nSeq = float(len(alignment)) # number of sequences in alignment
    consensus = '' # empty consensus string to be filled in and passed back
    
    for i in range(n): # loop through index of all positions (columns) within alignment
        counts = {} # each iteration we add the number of each type of residues observed at the position
        for seq in alignment: # iterating each position of each sequence
            letter = seq[i] # appropriate letter for current sequence given by the position index
            if letter == '-': # skip remainder of loop if there's a dash
                continue
            counts[letter] = counts.get(letter, 0) + 1 #increase count for the type of residue by1 via .get() so starting value of zero is obtained if not seen the letter before cycle
   
        fractions = []
        for letter in counts:
            frac = counts[letter]/nSeq # proportion of the total it represents
            fractions.append([frac, letter]) # list is used instead of dictionary so we can sort via numeric value as its first
    
        fractions.sort()
        bestFraction, bestLetter = fractions[-1] # want best to build consensus, last element = best as sorted low to high
    
        if bestFraction < threshold:
            consensus += 'X' # winning fraction below significance threshold extend with X
        else:
            consensus += bestLetter # could also put letters in a list and used ''.join()
    
    return consensus

alignment = ['SRPAPVVIILIILCVMAGVIGTILLISYGIRLLIK',
             'TVPAPVVIILIILCVMAGIIGTILLISYTIRRLIK',
             'HHFSEPEITLIIFGVMAGVIGTILLISYGIRRLIK',
             'HEFSELVIALIIFGVMAGVIGTILFISYGSRRLIK']
print(consensus(alignment))


###############################################################################
#####                       ALIGNMENT PROFILE                             #####
###############################################################################
# Profile: per position statistic of frequency of residue at given location
#   fractions of residue types, rather than single residues

def profile(alignment):
    
    n = len(alignment[0])
    nSeq = float(len(alignment))
    prof = []
    
    for i in range(n):
        counts = {}
        
        for seq in alignment:
            letter = seq[i]
            if letter == '-':
                continue
            counts[letter] = counts.get(letter, 0) + 1
            
        for letter in counts:
            counts[letter] /= nSeq
        prof.append(counts)
    
    return prof

alignment = ['SRPAPVVIILIILCVMAGVIGTILLISYGIRLLIK',
             'TVPAPVVIILIILCVMAGIIGTILLISYTIRRLIK',
             'HHFSEPEITLIIFGVMAGVIGTILLISYGIRRLIK',
             'HEFSELVIALIIFGVMAGVIGTILFISYGSRRLIK']
print(profile(alignment))
# output is a list of dictionaries, alignment positions are within a list index
# sub-dictionary for the index gives the fractions of residues present
# position-specific scoring matrices - like in PSI-BLAST
# POINT: build profile for family of sequences with commonality then searching
# other sequences with the whole profile, find seqeunces with shared properties
# of aligned family as a whole
# family profile displays presence of highly conserved or invariant site
# invariant sites would not be displayed if substitution table used
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
# PROFILE ALIGNMENTS: 
def profileAlign (profileA, profileB, simMatrix, insert=8, extend=4):
    
    numI = len(profileA) + 1
    numJ = len(profileB) + 1
    
    scoreMatrix = [[0] * numJ for x in range(numI)]
    routeMatrix = [[0] * numJ for x in range(numI)]
    
    for i in range(1, numI):
        routeMatrix[i][0] = 1
    for j in range(1, numJ):
        routeMatrix[0][j] = 2
    for i in range(1, numI):
        for j in range(1, numJ):
            
            penalty1 = insert
            penalty2 = insert
            
            if routeMatrix[i-1][j] == 1:
                penalty1 = extend
            elif routeMatrix[i][j-1] == 2:
                penalty2 = extend
            
            fractionsA = profileA[i-1]
            fractionsB = profileB[j-1]
            
            similarity = 0.0
            totalWeight = 0.0
            for residueA in fractionsA:
                for residueB in fractionsB:
                    weight = fractionsA[residueA] * fractionsB[residueB]
                    totalWeight += weight
                    similarity += weight * simMatrix[residueA][residueB]
            
            penalty1 *= totalWeight
            penalty2 *= totalWeight
            
            paths = [scoreMatrix[i-1][j-1] + similarity, # Route 0
                     scoreMatrix[i-1][j] - penalty1,     # Route 1
                     scoreMatrix[i][j-1] - penalty2]     # Route 2
            
            best = max(paths)
            route = paths.index(best)
            
            scoreMatrix[i][j] = best
            routeMatrix[i][j] = route
        
        profileOutA = []
        profileOutB = []
        
        i = numI-1
        j = numJ-1
        score = scoreMatrix[i][j]
        
        while i > 0 or j > 0:
            route = routeMatrix[i][j]
            
            if route == 0: # Diagonal
                profileOutA.append(profileA[i-1])
                profileOutB.append(profileB[j-1])
                i -= 1
                j -= 1
                
            elif route == 1: # Gap in profile B
                profileOutA.append(profileA[i-1])
                profileOutB.append(None)
                i -= 1
            
            elif route == 2: # Gap in profile A
                profileOutA.append(None)
                profileOutB.append(profileB[j-1])
                j -= 1
       
        profileOutA.reverse()
        profileOutB.reverse()
    
    return score, profileOutA, profileOutB

alignA = ['SRPAPVV--LII', 'TVPAPVVIILII']
alignB = ['HHFSEPEITLIIF', 'H-FSELVIALIIF']

profal = (profileAlign(profile(alignA), profile(alignB), BLOSUM62))
# output is numeric score and two aligned gapped profiles

# CONSENSUSMULTIPLEALIGN function given in documents for aligning based on 
# multiple consensus sequences

###############################################################################
#####                   PROFILE-BASED MULTIPLE ALIGNMENT                  #####
###############################################################################
# take list of unaligned sequences and a substitution matrix as input -
# initializing a value for the number of input sequences

from Alignments import BLOSUM62, sequenceAlign

def simpleProfileMultipleAlign(seqs, simMatrix):
    
    n = len(seqs)
    # align first two sequences, placing gapped output in multipleAlign list
    score, alignA, alignB = sequenceAlign(seqs[0], seqs[1], simMatrix)
    multipleAlign = [alignA, alignB]
    
    # loop through index for each sequence starting at 2 as we already used first
    # two sequences to get multiple alignment started - generate two new profiles
    # in loop, for existing alignment and for sequence to be added, seqs[i]
    for i in range(2, n):
        profA = profile(multipleAlign)
        toAdd = [seqs[i],]
        profB = profile(toAdd)
    
    # gives us alignment of two profiles and get two gapped profile alignments
    score, alignA, alignB = profileAlign(profA, profB, simMatrix)
    
    # use positions of gap insertions to combine new sequence into multiple alignment
    # output profiles are guides for placing insertions
    # repeat this operation twice, collecting and inserting gaps and once for newly added sequences
    gaps = []
    for j, fractions in enumerate(alignA): # enumerate loop for extraction of index number and fractions dictionary
        if fractions is None: # means we have a gap in profile
            gaps.append(j) # place current index in a list of gap positions
    # take the added column of '-' gap locations and put dashes into original multipleAlign list
    for j, seq in enumerate(multipleAlign):
        for gap in gaps:
            seq = seq[:gap] + '-' + seq[gap:]
            multipleAlign[j] = seq
    
    # second round gap insertion for locations from second align profile, placing gaps in toAdd list
    gaps = []
    for j, fractions in enumerate(alignB):
        if fractions is None:
            gaps.append(j)
    for j, seq in enumerate(toAdd):
        for gap in gaps:
            seq = seq[:gap] + '-' + seq[gap:]
        
        toAdd[j] = seq
    
    multipleAlign.extend(toAdd)
    
    return multipleAlign
    
# TEST
seqs = ['SRPAPVVLIILCVMAGVIGTILLISYGIRLLIK',
        'TVPAPVVIILIILCVMAGIIGTILLLIISYTIRRLIK',
        'HHFSEPEITLIIFGVMAGVIGTILLLIISYGIRLIK'
        'HFSELVIALIIFGVMAGVIGTILFISYGSRLIK']

align = simpleProfileMultipleAlign(seqs, BLOSUM62)
for k, seq in enumerate(align):
    print(k, seq)
    

###############################################################################
#####                       ClustalW from Python                          #####
###############################################################################
# Interfacing multiple alignments can be done with: ClustalW, MUSCLE, T-Coffee or BioPython modules

import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO, AlignIO

fastaFileName = "test2.fasta"
alignFileName = "testa.aln"

records = []
for i, seq in enumerate(seqs): # loop throughh sequences to make list of sequence record objects
    seqObj = Seq(seq, IUPAC.protein) # each one-letter sequence we make an object using standard IUPAC
    name = 'test%d' % i # combine IUPAC with name
    recordObj = SeqRecord(seqObj, id=name, description='demo only') # object for writing FASTA file
    records.append(recordObj)

outFileObj = open(fastaFileName, "w") # creates file handle object into which sequence records are written
SeqIO.write(records, outFileObj, "fasta") #writes in FASTA format
outFileObj.close()

# input file made, now want to run external alignment program - ClustalW
# would have to have ClustalW installed

from subprocess import call
cmdArgs = ['clustalw',
           '-INFILE=' + fastaFileName,
           '-OUTFILE=' + alignFileName]
call(cmdArgs)

# use .read() for one alignment in file, .parse() for several and getting back a list
fileObj = open(alignFileName)
alignment = AlignIO.read(fileObj, "clustal")
# reading functions makes alignment object from which we can print out attributes, and loop through to indiv records
print("Alignment length %i" % alignment.get_alignment_length())
for record in alignment:
    print(record.seq, record.id)

alignments = [alignment,]
outputHandle = open("test2.phylip", "w")
AlignIO.write(alignments, outputHandle, "phylip")
