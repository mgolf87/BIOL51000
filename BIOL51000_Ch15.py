###############################################################################
#####                 CHAPTER 15: MACROMOLECULAR STRUCTURE                #####
###############################################################################
# 3D arrangements of biological moleculues - 4D when we consider time-dependency and dynamics
# Chapter focues on - protein and directionaly functional untranslated RNA - fxn
# dependent on 3D structure
# Protein Structure: know atoms and covalent bonds, deviations are dynamic chemical
# groups like acidic residues (H on or off), enzyme modifications - post translational
# modifications - cross link cysteine, cut backbone, adding sugars, fats, phosphate groups
# Native conformation: generally lowest energy state and the range around it as molecule
# is not static due to temp and kinetic energy. Higher temp = larger range of native states
# covalent disulphide links form under oxidizing conditions, hydrophillic and phobic
# interactions establish compactness
# N-H amide and C=O carboxyl form polar hydrogen bonds
# hydrophillic exterior + hydrophobic core = globule
# amino acid sequences without hydrophobic components dont fold as no core and so
# are highly dynamic unstructured regions at the end of protein chains acting as
# flexible linkers between folded domains which are compact and globular

### PROTEIN STRUCTURES ###
# structural hierarchy as we can understand the final product as a sum of interactions
# or smaller elements
# PRIMARY: sequence of amino acids in the polypeptide chain including disulphide 
# links and post-translational modifications - covalent bond connectivity of the residues
# left handed biologically abundant chiral form present never right unless stated
# covalent bonds
# SECONDARY: the regular arrangment of hydrogen bonding along the backbone giving rise
# to the twist angles (dihedral/torsion) of the chain. Alpha-helix, beta-sheets, turn,
# and random coil (no regular structure). Ramachandran angles - twist of backbone relative
# to the alpha carbon - side chain branches
# hydrogen bonds
# TERTIARY: folded 3D shape and structure of one protein molecule in isolation
# QUATERNARY: structure resultant from multiple protein interactions forming a larger
# composite structure or complex (can include RNA, DNA, small molecules), can be same copies
# of proteins operating in highly symmetric structure

### MEMBRANE PROTEINS ###
# exclusion of water = differential cellular compartments = membrane proteins to
# transport compounds across the barrier
# transmembrane domains = generally inserted in membrane hole as they are made
# amino acid composition - special signal at start of sequence + lots of hydrophobic
# amino acids = transmembrane domain or a-helical bundles and B-barrels, 30% all
# proteins have a transmembrane domain - hard to isolate in experiments

### RNA ###
# hydrogen bonding of basepairs = complementary interactions = secondary structure
# fold back on itself, or interact with another RNA = duplex, stem-loops
# full conformation = tertiary = biological catalysts = ribozymes

### MACROMOLECULAR STRUCTURE ###
# variety of methods with varying precision and resolution for determining structure
# majority of high-resolution structures found with X-ray crystallography and nuclear
# magnetic resonance
# X-ray: form crystals of molecule of interest = most common structure determination
# growing crystals hardest part, molecules held in lattice = regular array of atom positions
# diffracts a beam of X-ray radiation = diffraction pattern = spots called reflections
# collect several diffraction patterns via vary angles = combine to determine the
# 3D map of electron density of the atom - fits covalent chemical structures to the density
# once best fit established it will do the same for the side chains
# NMR: second most common structure determination, smaller structures, size limit for study
# concentrated solution placed in electromagnet to generate strong magnetic field where
# the spin-active isotopes of the atomic nuclei will align to the field. H1 and P31 isotopes
# abundant in nature, sample enriched with C13 and N15 isotopes. Magnetic alignments detected
# by passing pulse of radio waves through the sample = creates resonance = frequency is dependent
# on the chemical and structural environment, moving pulses from atom to atom to determine which
# is present where = correlate observed resonances to connect together for 3D arrangment

### COMPARATIVE MODELING OR HOMOLOGY MODELING ###
# more data provided = better fit to native structure
# relies on the fact that when proteins evolve they change their aa seq more than their overall structure
# sequence similarity -> infer common ancestry -> structural similarities
# two basic steps for building protein homolog structure: find homologue of known
# structure and then use the homologues structure to guide building of a module
# query of unknown structure -> seq find homolog w/ known structure = structural template
# comparative modeling uses family-specific scoring matrices or substituion tables specific
# to the structural environment, structural environment for each position in sequence
# alignment is known based on the structural template. Environments defined by 
# combinig side chain H bonding, solvent exposure, secondary-structure categories
# examples: exposed alpha-helix, buried, side chain H bonds alpha helix, exposed beta sheets
# backbone built first, then side chains as vary, then loops modelled in unaligned regions = gaps
# Alternative: spatial restraints derived from the templates on the model of the query
# minimisattion performed after to find conformation to best satisfi the restraints = MODELLER


###############################################################################
#####                     OBTAINING STRUCTURE DATA                        #####
###############################################################################
# power of Python to download data directly from PDB websites download service
from urllib.request import urlopen
# URL for PDB data where %s identifying database
PDB_URL = 'http://www.rcsb.org/pdb/cgi/export.cgi/' \
            '%s.pdb?format=PDB&compression=None'
def downloadPDB(pdbId, fileName=None):
    if not fileName:
        fileName = '%s.pdb' % pdbId
        
    response = urlopen(PDB_URL % pdbId)
    data = response.read().decode('utf-8') # converts to string via decoding as imported as bytes
    fileObj = open(fileName, 'w')
    fileObj.write(data)
    fileObj.close()
    
    return fileName
fileName = downloadPDB('1A12')

###############################################################################
#####                    GEOMETRIC MANIPULATION                           #####
###############################################################################
# First function to find the center of mass or centroid of a structure
from numpy import zeros
# dictionary of the atomic numbers of common atoms in biological materials
ATOMIC_NUMS = {'H':1, 'C':12, 'N':14, 'O':16, 'P':31, 'S':32}
# increase accuracy with more elements and masses rather than numbers

def getCenterOfMass (structure):
    centerOfMass = zeros(3, float) # initialize numpy.array object of zeros
    # 3 coordinate positions to define the center in 3D
    totalMass = 0.0 # set 0 to initialize
    
    for chain in structure.chains: # iterate over all chains
        for residue in chain.residues: # iterate over all residues
            for atom in residue.atoms: # iterate over all atoms
                mass = ATOMIC_NUMS.get(atom.element, 12.0) # atom.element as key to dict, multiply to further differentiate from H and other element masses
                centerOfMass += mass * atom.coords
                totalMass += mass
    
    centerOfMass /= totalMass # per atom average
    
    return centerOfMass
# TEST #
struc = getStructuresFromFile('examples/1A12.pdb')[0]
print(getCenterOfMass(struc))

### CHANGE THE STRUCTURAL COORDINATES ###
# translate and rotate a structure without fucking the initial coords - rotation matrix
# 3D space a rotation matrix is 3x3 square array of numbers representing transformation
# of Cartesian coordinates
from numpy import array, dot
from Maths import getRotationMatrix
from math import pi

def rotateStructure(structure, axis=(1,0,0), angle=0):
    rMatrix = array(getRotationMatrix(axis, angle))
    
    for chain in structure.chains:
        for residue in chain.residues:
            for atom in residue.atoms:
                newCoords = dot(rMatrix, atom.coords)
                atom.coords = newCoords
rotateStructure(struc, (0,1,0), pi/2)
# each atom calc dot produce of rotation matrix and original coordinate array

###############################################################################
#####                   STRUCTURES --> NUMPY ARRAYS                       #####
###############################################################################

from numpy import array, dot
def getAtomCoords (structure):
    
    # initialize new numpy objects
    coords = []
    atoms = []
    
    # loop through all chains, residues, and atoms in the stucture adding to 
    # the initialized numpy objects. We know which coordinate for each atom
    # hierarchical molecular data -> numeric numpy array
    for chain in structure.chains:
        for residue in chain.residues:
            for atom in residue.atoms:
                coords.append(atom.coords)
                atoms.append(atom)
    
    return atoms, array(coords)

###############################################################################
#####                        DISTANCES AND ANGLES                         #####
###############################################################################



###############################################################################
#####                          STRUCTURAL SUBSETS                         #####
###############################################################################

# weights paramter allows for different atoms to have different degrees of influence
# generally due to their relative atomic mass
def centerCoords(coords, weights):
   
    # transpose is done to ensure multiplication with three rows separately with weights
    wCoords = coords.transpose() * weights
    # sum along the rows [axis=1] -> array of three numbers
    xyzTotals = wCoords.sum(axis=1) 
    center = xyzTotals/sum(weights) # produces an average or center of coordinates
    coords -= center # move to new center by removing old center position
    
    return coords, center

###############################################################################
#####                            COORDINATES                              #####
###############################################################################
import numpy as np
from numpy import zeros, ones, cross, sqrt, linalg, exp, identity
def alignCoords(coordsA, coordsB, weights=None):
    """ Takes two arrays of coordinates, optional array of weights, and finds
        the rotation which best superimposes pairs of coordinates """
    
    n = len(coordsA)
    if weights is None: # generates series of 1's same length as coordinates
        weights = ones(n)
        
    # find an optimum rotation transformation to minimize the differences in positions
    # between corresponding coordinates    
    rMat = dot(coordsB.transpose()*weights, coordsA) # weighted dot product of one 
    # coordinate array and the transpose of another = transformation of positions
    
    # Singular value decomposition (SVD) - factorises matrix into components
    # components are a rotation matrix, array of linear scaling factors, and an
    # opposing rotation matrix
    rMat1, sclaes, rMat2 = linalg.svd(rMat)
    
    # need to determine if mirror image so we multiply the detriments
    sign = linalg.det(rMat1) * linalg.det(rMat2)
    # detriment product determines our next action - if mirror image or not
    if sign < 0:
        rMat1[:, 2] *= -1 # flips the sign = remove reflection
    rotation = dot(rMat1, rMat2) # optimised rotation matrix
    coordsB = dot(coordsB, rotation) # apply coordinate transformation
    
    return rotation, coordsB


###############################################################################
#####                     ROOT MEAN SQUARE DEVIATION                      #####
###############################################################################
# Calculate variation in coordinates across the atom positions represented in arrays
# differences in coordinate positions, square them, find average and square the root
# average distance of coordinate spread -> average of squares not distances 
# = bias towards larger deviations having more influence
from math import sqrt 

def calcRmsds (refCoords, allCoords, weights):
    """ takes array of reference coordinates, list of the other coordinate
    array for comparison and a list of weights to bias atoms separately """
    rmsds = [] # hold values
    totalWeight = sum(weights) # total of all input weights
    totalSquares = zeros(refCoords.shape) # empty array of same size to hold summation of positional differences
    
    # loop through coordinate arrays comparing to reference
    # whole array object manipulations - applied to all elements of the array
    # delta = array of all coordinate differences, squaring the elements
    for coords in allCoords:
        delta = coords - refCoords
        squares = delta * delta
        totalSquares += squares
        sumSquares = weights*squares.sum(axis=1) # squared deviation for each atom
        # sum of square values along the spatial axis. multiply by weights for sumSquares
        rmsds.append(sqrt(sum(sumSquares)/totalWeight))
    
    nStruct = len(allCoords) 
    atomRmsds = sqrt(totalSquares.sum(axis=1)/nStruct) # average of the RMSDs
    # for each atoms RMSD over the whole set of conformations
    
    return rmsds, atomRmsds
# finds optimum rotation to superimpose two arrays of coordinates

###############################################################################
#####                     ENSEMBLE STRUCTURE ALIGNMENT                    #####
###############################################################################
# We want to superimpose more than two coordinate arrays -> superimpose all conformations
# of a whole structure ensemble via repeating pairwise superimpositions relative to reference.
# A better method would compare all against all at the same time - quantum programming?
def superimposeCoords(allCoords, weights, threshold=5.0):
    nStruct = len(allCoords)
    
    refCoords = allCoords[0] # set to first element of list of arrays
    meanCoords = zeros(refCoords.shape) # used to add coordinates after aligned to reference
    rotations = [] # starts empty and is passed back at the end once filled with alignments
    
    # first superimposition = loop all coordinate arrays aligning to reference
    for index, coords in enumerate(allCoords):
        if index == 0: # don't need to align with self as its the reference
            rotation = identity(3) # no rotation
        else: # coordinate alignment
        # optimised rotation matrix is updated and put back into allCoords
            rotation, coords = alignCoords(refCoords, coords, weights)
            allCoords[index] = coords # update to aligned - replaces previous
        
        rotations.append(rotation) # rotation matrix collected into rotations list
        meanCoords += coords # coordinate array added to total
    
    meanCoords /= nStruct # gives us the average positions
    rmsds, atomRmsds = calcRmsds(meanCoords, allCoords, weights) # calcs RMSDs and stores
    
    bestRmsd = min(rmsds) # best - or minimum value (smallest)
    bestIndex = rmsds.index(bestRmsd) # closest to mean - min index
    bestCoords = allCoords[bestIndex] # best set of coords now reference set for another round
    # bestCoords = reference array closest to mean of aligned coordinates
    
    # adjusted weights for second round scaled by threshold value
    # new weights have atoms with small RMSD values having largest weights
    # and as the variance increases the weighting reduces exponentially
    weightScale = atomRmsds/threshold
    weights *= exp(-weightScale*weightScale)
    
    meanCoords = bestCoords.copy() # average coordinates = closest to mean
    
    # skip indexes which match reference coords, and for all other indexes
    # superimposition alignment executed = optimising rotation to reference
    # we are inserted updated coords into allCoords at their index to then find
    # the new coord average after completion
    for index, coords in enumerate(allCoords):
        if index != bestIndex:
            rotation, coords = alignCoords(bestCoords, coords, weights)
            rotations[index] = rotation
            allCoords[index] = coords # updates to aligned
    
    meanCoords /= nStruct # new average coord set
    
    # compare new coords with average set of coords
    # new unbiased weights set to 1 via 'ones' = observed distance variation for each
    # atom is reported equally
    weights = ones(len(weights))
    rmsds, atomRmsds = calcRmsds(meanCoords, allCoords, weights)
    
    # optimized superimposition of coordinate arrays done
    return allCoords, rmsds, atomRmsds, rotations

# applies the optimised superimposition of coordinate arrays to structure objects
def superimposeStructures(structures):
    """ Takes a list of structure objects extracting coordinate arrays, 
        performs coordinate superimposition, and updates original Atom objects
        with new coordinates """
        
    weights = None # initialize
    allCoords = [] # all coordinates for all structures - list of 2D arrays
    
    # structure inspection - generating coordinate arrays and atom object lists
    for structure in structures:
        atoms, coords = getAtomCoords(structure)
        # defines weight list with atomic numbers - not very precise but fine
        if weights is None:
            weights = [ATOMIC_NUMS[atom.element] for atom in atoms]
            weights = array(weights)
        
        # redefine coordinates via moving to center
        coords, center = centerCoords(coords, weights)
        allCoords.append(coords) # centered array put into list containing all coords
        # run function to superimpose
        results = superimposeCoords(allCoords, weights)
        # Results = modified array coordinates, list of RMSD values for each 
        # structure and atom, list of rotations for the superimposition
        allCoords, rmsds, atomRmsds, rotations = results
    
    # generate list index and structure object from the structures
    for i, structure in enumerate(structures):
        atoms, oldCoords = getAtomCoords(structure) # list of atoms and array of original coords
        
        # update atom coords for each atom object to the new superimposed coords
        for j, atom in enumerate(atoms):
            atom.coords = allCoords[i][j]
            
    return rmsds, atomRmsds

###############################################################################
#####                 HOMOLOGOUS STRUCTURE ALIGNMENT                      #####
###############################################################################



###############################################################################
#####                               BIOPYTHON                             #####
###############################################################################
# load PDB module then create parser object to make PDB.Structure object using data from file
from Bio import PDB
fileName = 'examples/1UST.pdb'
parser = PDB.PDBParser()
struc = parser.get_structure('Name', fileName)

# extract first set of coordinates (first conformation) and loop through all chains
# residues and atoms to get coordinates
conformation = struc[0]
for chain in conformation:
    for residue in chain:
        atomNames = [a.name for a in residue]
        print(chain.id, residue.id[1], residue.resname, atomNames)
        caAtom = residue['CA']
        print(caAtom.name, caAtom.coord, caAtom.bfactor)
# writing data to disk via writer
outFileName = 'test.pdb'
writer = PDB.PDBIO()
writer.set_structure(struc)
writer.save(outFileName)

def findCloseAtomsBioPy (structure, xyz, conf=0, limit=5.0):
    
    if not isinstance(structure, PDB.Structure.Structure):
        raise Exception('Structure must be Bio.PDB.Structure class')
    
    closeAtoms = []
    xyz = array(xyz)
    limit2 = limit * limit
    
    coords = []
    atoms = []
    confModel = structure[conf]
    
    for chain in confModel:
        for residue in chain:
            for atom in residue:
                coords.append(atom.coord)
                atoms.append(atom)
    deltas = coords - xyz
    squares = deltas * deltas
    sumSquares = squares.sum(axis=1)
    
    boolArray = sumSquares < limit2
    indices = boolArray.nonzero()[0]
    closeAtoms = [atoms[i] for i in indices]
    
    return closeAtoms

###############################################################################
#####                              PYMOL                                  #####
###############################################################################
import pymol
pymol.finish_launching() # graphical environment to appear
# PDB file previously downloaded can then be loaded
fileName = downloadPDB('examples/1AFO')
strucName = 'Glycophorin'
pymol.cmd.load(fileName, strucName)
# define two subsets of structure 'bb' (backbone) and 'tmd' (transmembrane)
pymol.cmd.select('bb', 'name n+c+o+ca')
pymol.cmd.select('tmd', 'resi 71-100 and name n+c+o+ca ')
# issue commands for structure display
pymol.cmd.color('gray', strucName)
pymol.cmd.color('red', 'bb')
pymol.cmd.show('cartoon', 'bb')
pymol.cmd.color('blue', 'tmd')
pymol.cmd.hide('lines', 'hydro')

outFileName = strucName + '.pdb'
pymol.cmd.save(outFileName, strucName, 0, 'pdb')
pymol.cmd.png(strucName+'.png')
pymol.cmd.quit()