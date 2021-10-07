###############################################################################
#####                      CHAPTER 16: ARRAY DAY                          #####
###############################################################################
# represent MA data as 2D array in numpy

from numpy import array, dot, log, sqrt, uint8, zeros

### Loading data from file and creating Microarray object
###############################################################################
##################### Importing Text Matrices #################################
###############################################################################
def loadDataMatrix(fileName, sampleName, default=0.0):
    
    fileObj = open(fileName, 'r')
    
    """ Empty sets for identifiers of row and cols in microarray data
    identifiers = array coordinates or text labels -> hashable
    Values will be keys to access numeric dataDict data """
    rows = set()
    cols = set()
    dataDict = {}
    
    for line in fileObj:
        row, col, value = line.split() # split on whitespace
        if row not in dataDict: # checks for row identifiers having entry
            dataDict[row] = {} # make new inner dictionary with row key
        if col in dataDict[row]:
            print('Repeat entry found for element %d, %d' % (row, col))
            continue
        
        dataDict[row][col] = float(value) # actual signal value w/ r,c key
        rows.add(rows) # added to the set
        cols.add(col) # sets ignore repeat cases
        
        rows = sorted(rows) # sort to ordered lists
        cols = sorted(cols) # total range of identifiers collected
        
        nRows = len(rows) # used to create axes of 2D array
        nCols = len(cols)
        
        dataMatrix = zeros((nRows, nCols), float) # numpy array of zeros
        # fill array via extracting values from dataDict - replace missing w/ default
        for i, row in enumerate(rows): # extract indexes as we loop identifiers
            for j, col in enumerate(cols):
                value = dataDict[row].get(col, default)
                dataMatrix[i,j] = value # filling based on extracted indexes
        
        fileObj.close()
    
    return Microarray(sampleName, dataMatrix, rows, cols)

### Assignment 4 - 1/4 - change 3 to 2>?
###############################################################################
##################### Extracting Array Image Data #############################
###############################################################################
from PIL import Image
from Images import imageToPixmapRGB
from numpy import array, dstack, transpose, uint8, zeros, log2

def loadArrayImage (fileName, sampleName, nRows, nCols=None):
    
    """ Take signal from each spot as total amount of signal within
    each grid cell """
    
    if not nCols: # nCols not given its set to nRows
        nCols = nRows
    
    # array of 0's with 3 layers to store color components
    dataMatrix = zeros((3, nRows, nCols), float) 
    
    img = Image.open(fileName) # automatic file type
    pixmap = imageToPixmapRGB(img) # converts to a numeric array
    
    height, width, depth = pixmap.shape
    
    dx = width/float(nCols) # used for precise values for grid start points (floats)
    dy = height/float(nRows) # These give us grid size
    xSize = 1 + (width-1)//nCols # (int) for for end points = fixed # pixels
    ySize = 1 + (height-1)//nRows # +1 pixel for slice of array, -1 so we dont overshoot edge
    
    for row in range(nRows): 
        yStart = int(row*dy) # first pixel position for image section (row#*rowDepth)
        yEnd = yStart + ySize # last pixel be inside the limit
        for col in range(nCols): # within each row -> calc range of pixels to
            xStart = int(col*dx) # select a column of data from the image
            xEnd = xStart + xSize
    
    elementData = pixmap[yStart:yEnd, xStart:xEnd]
    dataMatrix[:, row, col] = elementData.sum(axis=(0,1)) # gives total signal for grid element 
    # : is for the first axis -> setting al color channels at the same time
    
    return Microarray(sampleName, dataMatrix)

###############################################################################
#####                       MICROARRAY CLASS                              #####
###############################################################################
class Microarray(object):
    
    # called each time object created -> (Sample Name, NumPy Array)
    # can also give names for row and col data labels
    def __init__(self, name, data, rowData=None, colData=None):
        
        self.name = name # stores name for MA data, self -> object called
        data = array(data) # copy original and convert lists or tuples
        
        shape = data.shape # sizes of axes extracted
        
        if len(shape) == 3: # represent colours
            self.nChannels, self.nRows, self.nCols = shape
        
        elif len(shape) == 2: # one channel
            self.nRows, self.nCols = shape
            self.nChannels = 1 # forcing 3 dimensions 
            data = array([data])
        
        else: # array data axes dont have 2 or 3 axes
            raise Exception("Array data must have either 2 or 3 axes.")
        
        self.data = data # adjusted data saved to object
        self.origData = array(data) # original numpy data
        
        self.rowData = rowData or range(self.nRows) # labels associated w/ object
        self.colData = colData or range(self.nCols) # either list of labels is None
                                                    # will be defined as range
    # Constructor Function End #
    def resetData(self): # enables us to revert to original data if needed
        self.data = array(self.origData) 
        self.nChannels = len(self.data)
    
    # Special Functionality #
    def writeData(self, fileName, separator=' '):
        # opens file to write out
        fileObj = open(fileName, 'w')
        
        # loop through rows and cols converting identifiers
        # to strings from numbers
        for i in range(self.nRows):
            rowName = str(self.rowData[i])
            
            for j in range(self.nCols):
                colName = str(self.colData[i])
                
                values = self.data[:, i, j] # (i,j) to get data from array for all array channels
                
                # text written to file -> list of data of name of row/col at start
                lineData = [rowName, colName] # w/ string representations of numeric
                lineData += ['%.3f' % (v,) for v in values] # convert floats to strings w/ 3 decimals
                
                line = separator.join(lineData) # combines separate lineData strings
                fileObj.write(line + '\n') # write the line to file object
    
    def makeImage (self, squareSize=20, channels=None):
        
        minVal = self.data.min()
        maxVal = self.data.max()
        dataRange = maxVal - minVal
        
        adjData = (self.data - minVal) * 255 / dataRange # scales the colors
        adjData = array(adjData, uint8)
        
        if not channels:
            if self.nChannels == 1:
                channels = (0,0,0) # Greyscale
            else:
                channels = list(range(self.nChannels))[:3]
            
        pixmap = []
        for i in channels:
            if i is None:
                pixmap.append(zeros((self.nRows, self.nCols), uint8))
            else:
                pixmap.append(adjData[i])
                
        while len(pixmap) < 3:
            pixmap.append(zeros((self.nRows, self.nCols), uint8))
        
        pixmap = dstack(pixmap) # stacks along depth axis
        img = Image.fromarray(pixmap, 'RGB')
        
        width = self.nCols * squareSize
        height = self.nRows * squareSize
        img = img.resize((width, height))
        
        return img
    
    def clipBaseline (self, threshold=None, channels=None, defaultProp=0.2):
        
        if not channels:
            channels = range(self.nChannels)
        
        channels = [tuple(channels)]
        
        maxVal = self.data[channels].max()
        if threshold is None:
            limit = maxVal * defaultProp
        else:
            limit = threshold
        
        boolArray = self.data[channels] < limit
        indices = boolArray.nonzero()
        
        self.data[indices] = limit
        
        self.data[channels] -= limit
        self.data[channels] *= maxVal / (maxVal-limit)
    
    def normaliseSd (self, scale=1.0):
        for i in range(self.nChannels):
            self.data[i] = self.data[i] * scale / self.data[i].std()
    
    def normaliseMean(self, scale=1.0):
        for i in range(self.nChannels):
            self.data[i] = self.data[i] * scale / self.data[i].mean()
    
    def centerMean(self):
        for i in range(self.nChannels):
            self.data[i] -= self.data[i].mean()
    
    def normaliseZscore(self):
        self.centerMean()
        self.normaliseSd()
        
    def normaliseMax(self, scale=1.0, perChannel=True):
        if perChannel:
            for i in range(self.nChannels):
                self.data[i] = self.data[i] * scale / self.data[i].max()
        else:
            self.data = self.data * scale / self.data.max()
    
    def normaliseRowMax(self, scale=1.0):
        for i in range(self.nChannels):
            self.data[i] = self.data[i] * scale / self.data[i].max(axis=1) [:, None]
    
    def normaliseRowMean(self, scale=1.0):
        for i in range(self.nChannels):
            self.data[i] = self.data[i] * scale / self.data[i].mean(axis=1) [:, None]
    
    def normaliseColMax (self, scale=1.0):
        for i in range(self.nChannels):
            self.data[i] = self.data[i] * scale / self.data[i].max(axis=0)
    
    def normaliseColMean(self, scale=1.0):
        for i in range(self.nChannels):
            self.data[i] = self.data[i] * scale / self.data[i].mean(axis=0)
    
    def normaliseRefs (self, rows, cols, scale=1.0, channels=None):
        if not channels:
            channels = range(self.nChannels)
        channels = tuple(channels)
        refValues = self.data[channels, rows, cols]
        
        for i in channels:
            self.data[i] = self.data[i] * scale / refValues[i].mean()
    
    ### Normalize fluorescence intensity via log of mean
    def normaliseLogMean(self):
        self.clipBaseline(threshold=0.0)
        for i in range(self.nChannels):
            self.data[i] = log(1.0 + self.data[i] / self.data[i].mean())
            
    ###
    
    def normaliseQuantile(self, refData, channel=0):
        # could be to different channel
        
        values = self.data[channel].flatten()
        order = values.argsort()
        refValues = refData.flatten()
        refValues.sort()
        
        refSelection = order.argsort()
        values = refValues[refSelection]
        
        self.data[channel] = values.reshape((self.nRows, self.nCols))
    
    def normaliseRowQuantile(self, channel=0):
        
        channelData = self.data[channel]
        orders = channelData.argsort(axis=1)
        sortedRows = array(channelData)
        sortedRows.sort(axis=1)
        refValues = sortedRows.mean(axis=0) # average over columns
        
        rows = range(self.nRows)
        self.data[channel,rows,:] = refValues[orders[rows,:].argsort()]
        
    ### Changing Array Channels ###
    def checkDataSize(self, channelData):
        
        channelData = array(channelData)
        if channelData.shape != (self.nRows, self.nCols):
            msg = 'Attempt to use data of wrong size'
            raise Exception(msg)
        return channelData
    
    def setChannel(self, channelData, index=0):
        channelData = self.checkDataSize(channelData)
        self.data[index] = channelData
   
    def addChannel(self, channelData):
        from numpy import append
        channelData = self.checkDataSize(channelData)
        
        self.data = append(self.data, channelData, axis=0)
        self.nChannels += 1
    
    def combineChannels(self, indexA, indexB, combFunc=None, replace=None):
        if not combFunc:
            import operator
            combFunc = operator.add
        
        channelData = combFunc(self.data[indexA], self.data[indexB])
        
        if replace is None:
            self.addChannel(channelData)
        else:
            self.setChannel(channelData, replace)
    
    
###############################################################################
#####                 DIFFERENCES AND SIMILARITIES                        #####
###############################################################################
### Assignment 4 - Part 1
# Differences - G-Scores, gives greyscale
############################################ All part of MA Class
    imgFile = 'RedGreenArray.png'
    rgArray = loadArrayImage(imgFile, 'TwoChannel', 18, 17)
    diff = rgArray.data[0]-rgArray.data[1] # comparing two signal intensity arrays
# store the differences in the first two color channels
# flip sign for green channel and store, remove any negative
    rgArray.setChannel(diff, 0)
    rgArray.setChannel(-diff, 1)
    rgArray.clipBaseline(threshold=0.0, channels=(0,1))
    rgArray.makeImage(20).show()

# Similarities - if R/G gives Yellow
    from operator import mul # multiply

    rgArray = loadArrayImage(imgFile, 'TwoChannel', 18, 17)
    rgArray.combineChannels(0, 1, combFunc=mul, replace=2)
    rgArray.makeImage(20, channels=(2,2,None)).show()

### Assignment 4 - part 2
    #from numpy import log2
# Takes two input data arrays, produces combined comparison array
    def log2Ratio(data1, data2):
        data1 = array(data1) + 1e-3
        data2 = array(data2) + 1e-3
    
        return log2(data1/data2)

    rgArray = loadArrayImage(imgFile, 'TwoChannel', 18, 17)
    rgArray.combineChannels(0, 1, combFunc=log2Ratio, replace=2)

# G-Score -- log scaled
### Assignment 4 - part 4
# differences b/w two channel intensities via logs: x*log(x/y) where x and y are two color channels
# shows information content of one distribution over another - G-score
    #from numpy import log2
    def gScore(data1, data2):
        data1 = array(data1) + 1e-3
        data2 = array(data2) + 1e-3
        
        return data1 * log2(data1/data2)
# normalise values so logs are scaled into same range as other channels
    rgArray = loadArrayImage(imgFile, 'TwoChannel', 18, 17)
    rgArray.combineChannels(0, 1, combFunc=gScore, replace=2)
    rgArray.normaliseMax(perChannel=True)

    rgArray.makeImage(20, channels=(2, 2, 2)).show()


###############################################################################
#####                    HIERARCHIAL CLUSTERING                           #####
###############################################################################
# rearranging data (swapping rows and columns) to see similarities or correlations in elements
# same concept as phylo tree = most similar row/col placed into sub-groups

    def __hierarchialRowCluster(self, dataMatrix):
        """ Takes an array of data as 2D matrix to build distance matrix
        distance = similarity b/w rows
        measure = EucDist b/w row vectors (root of sum of differences squared)
        the __ creates a private function to be used internally w/in class"""
    
        from SeqVariation import neighbourJoinTree
    
        n = len(dataMatrix[0])
    
        distanceMatrix = zeros((n, n), float) # distance b/w row i and row j
    
        for channelData in dataMatrix: # loop through each layer and row
            for i, row in enumerate(channelData):
                diffs = channelData - row # subtract row from whole array
                sqDiffs = diffs * diffs # square the difference
                sqDists = sqDiffs.sum(axis=1) # sum along the row
                distanceMatrix[i,:] += sqDists # dist from row to all other rows
    
    # Creates hierarchial tree from distanceMatrix - list of lists
        tree, joinOrder = neighbourJoinTree(distanceMatrix.tolist())
    
        rowOrder = list(tree) # copy of the tree - list of nodes = row order
    
        i = 0
        while i < len(rowOrder): # size of rowOrder (sub-list) grows as branches are flattened
            while not isinstance(rowOrder[i], int): # checks if [i] is an int (=end of branch)
                                                    # and that its not in the sub-list rowOrder
                rowOrder[i:i+1] = rowOrder[i] # replacing original range covering the sub-list
        # second while loop needed -> new element at position i might be a sub-list
                i += 1
    
        return rowOrder


    def hierarchialCluster(self):
        # not private, clusters rows in self.data -> reorders via row hierarchy
        # transpose result (flip rows/cols) cluster cols (new rows) forming new col order
        rows = self.__hierarchicalRowCluster(self.data)
    
        swapped = transpose(self.data, axes=(0,2,1))
        cols = self.__hierarchicalRowCluster(swapped)
    
        data = self.data[:,rows] # rearrange
        data = data[:,:,cols]
    
        # data = array(data.tolist()) # to fix PIL.Image bug
    
        name = self.name + '-Sorted'
        rowData = [self.rowData[i] for i in rows]
        colData = [self.colData[j] for j in cols]
    
        sortedArray = Microarray(name, data, rowData, colData)
    
        return sortedArray






