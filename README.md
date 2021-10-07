# BIOL51000 - Data Systems in the Life Sciences
Lewis University Fall 2, 2020

<br />

### MGolf_BIOL51000_A3.ipynb 
Week 3 Assignment: Measuring Conservation, Calculating Substitution Matrices, Calculating Distance Matrices
### MGolf_BIOL51000_A4.ipynb 
Week 4 Assignment: Microarray Analysis
### MGolf_BIOL51000_A5.ipynb 
Week 5 Assignment: Calculating RMSD and Aligning a Structure Ensemble
### MGolf_BIOL51000_A6.ipynb 
Week 6 Assignment: Image Processing and GUIs

<br />

## Python for Biology, Tim Stevens
### BIOL51000_Ch12.py - Pairwise Sequence Alignments
Sequence Identity; Substitutability; Similarity; Pairwise Alignment; BLAST
### BIOL51000_Ch13.py - Multiple Sequence Alignments
Consensus; Alignment Profile; Profile-Based Multiple Alignment; ClustalW from Python
### BIOL51000_Ch15.py - Macromolecular Structure
Obtaining Structure Data; Geometric Manipulation; Structures -> Numpy Arrays; Structural Subsets; Coordinates; Root Mean Square Deviation; Ensemble Structure Alignment; BioPython; PYMOL
### BIOL51000_Ch16.py - Microarray Data
Importing Text Matrices; Extracting Array Image Data; Microarray Class; Differences and Similarities; Hierarchial Clustering

<br />

### Alignments.py
calcSeqIdentity; calcSeqSimilarity; pairAlignScore; sequenceAlign; makeBlastDatabase; blastSearch; makeBlastDatabase
### MultipleAlign.py
consensus; profile; profileAlign; simpleProfileMultipleAlign; consensusMultipleAlign
### Clustering.py
euclideanDist; findNeighbours; simpleCluster; dbScanCluster; kMeans; kMeansSpread; jumpMethodCluster; principleComponentAnalysis; extractPrincipleComponent; twoClassLda
### Files.py
readFastaFile (fragments); readFastaFile (sequence); calcCentroid; printPubmedAbstracts; writeSingleFastaSequence; writeFastaSeqs; writeListFile; writeChromosomeRegions; readListFile; readChromosomeRegions; writeCsvFile; readCsvFile; findFiles; removeFiles
### HTSequences.py
downloadFile; downloadGenomeFiles; uncompressGzFile; indexGenome; genomeAlign
### MachineLearning.py
getFeatureDistance; kNearestNeighbour; kNearestNeighbour (numpy); selfOrganisingMap; neuralNetPredict; neuralNetTrain; convertSeqToVector; kernelGauss; kernelLinear; svmTrain; svmPredict; svmSeparation
### Probability.py
binomialProbability; getNextGenPop; Viterbi; forwardBackward
### Sequences.py
proteinTranslation; estimateMolMass; matchDnaProfile; calcGcContent; hydrophobicitySearch; relativeEntropySearch; estimateCharge; estimateIsoelectric; calcRelativeEntropy
### SeqVariation.py
getConservation; makeSimilarityString; getAlignProperties; calcSubstitutionMatrix; getDistanceMatrix; getJoinPair; getDistToJunction; neighbourJoinTree; treeProfileMultipleAlign; ancestorResidue; calcSubstitutionRates; calcActivePassive
### Statistics.py
getMedian; binomialTailTest; poissonTailTest; normalTailTest; zTestMean; tConfInterval
### Structures.py
downloadPDB; writeStructureToFile; getCenterOfMass; rotateStructure; getAtomCoords; rotateStructureNumPy; affineTransformStructure; findCloseAtoms; getPhiPsi; filterSubStructure; copyStructure; alignCoords; calcRmsds; superimposeCoords; superimposeStructures; getChainSequence; seqStructureBackboneAlign; findCloseAtomsBioPy

<br />

#### .py files are adapted from Python for Biology, Tim Stevens
