""" Quad-tree for FMM
Observation: bit level operation can be optimized.
String variables not required.
"""

from numpy import *

############################################################
## Low level utilities to generate and manipulate bit numbers
############################################################
    
def toBinStr(n, level):
    ''' Convert a float to a binary string'''
    bStr=''
    if n < 0 or n > 1 : raise ValueError, 'number must be between [0,1]'
    if level < 0: raise ValueError, 'level must be a positive integer'
    levelcount = 0
    boxsize = 1.0         
    while levelcount < level:
        boxsize = boxsize / 2
        if n > boxsize:
            bStr = bStr + '1'
            n = n - boxsize
        else:
            bStr = bStr + '0'
        levelcount = levelcount + 1
    return bStr
    
def toBin(n, level):
    '''convert a float number to binary number... if printed it d be
    seen as a integer number'''
    bStr = toBinStr(n, level)
    return int(bStr, 2)
    
def binToStr(n):
    '''convert a binary number to a binary string'''
    bStr = ''
    while n > 0:
        bStr = str(n % 2) + bStr
        n = n >> 1
    return bStr
    
def bitInterleavingStr(n1, n2):
    '''apply bit interleaving to two binary numbers returning a binary 
    string'''
    bStr = ''
    if n1 + n2 == 0:
        # Especial case to manage the interleaving of two numbers zero
        bStr = '0'
    else:
        # General case to manage the interleaving of any positive number
        while n1 > 0 or n2 > 0:
            bStr = str(n1 % 2) + str(n2 % 2) + bStr
            n1 = n1 >> 1
            n2 = n2 >> 1
    return bStr
    
def bitInterleaving(n1, n2):
    '''apply bit interleaving to two binary numbers returning a binary 
    number'''
    bStr = bitInterleavingStr(n1, n2)
    return int(bStr, 2)
    
def bitDeinterleavingStr(n, l):
    '''apply bit deinterleaving to a number up to level l, returning the 
    two original binary numbers'''
    bStr1 = ''
    bStr2 = ''
    while n > 0 or l > 0:
        bStr2 = str(n % 2) + bStr2
        n = n >> 1
        bStr1 = str(n % 2) + bStr1
        n = n >> 1
        l = l - 1
    return bStr1, bStr2

def bitDeinterleaving(n, l):
    '''apply bit deinterleaving to a number, returning the two original
    binary numbers'''
    bStr1,bStr2 = bitDeinterleavingStr(n, l)
    return int(bStr1, 2), int(bStr2, 2)
    
############################################################
## Give high level access to Father, childs and Neighbors 
############################################################

def findFather(n):
    '''find the father s number of the box n'''
    return n >> 2
    
def findChilds(n):
    '''find the child s number of the box n'''
    shiftN = n << 2
    return [shiftN + 0, shiftN + 1, shiftN + 2, shiftN + 3]
    
def findNeighbors(n, l):
    '''find the neighbors numbers of the box n at level l'''
    # # Deinterleaving the number
    n1, n2 = bitDeinterleaving(n, l)
    
    # # Interleaving
    # create sets of neighbors, taking care of the border of the boxes
    # Creation of neighbors list combining the sets
    neighborSet = []
    for x in range(max(n1 - 1, 0), min(n1 + 1, 2**l - 1) + 1):
        for y in range(max(n2 - 1, 0), min(n2 + 1, 2**l - 1) + 1):
            if not (x == n1 and y == n2):
                neighborSet.append(bitInterleaving(x, y))

    return neighborSet
    
def findParentsChildren(n, l):
    '''
    Find the node s interaction list
    '''
    import operator
    return reduce(operator.add, [findChilds(pN) for pN in findNeighbors(findFather(n), l-1)])
    
############################################################
## Enable the use of objects as nodes of a quadtree
############################################################

class QuadtreeNode:
    '''
    Node of a quadtree for storing information
    '''
    
    def __init__(self, numberArg, levelArg):
        '''initialize the node'''
        self.level = levelArg
        self.number = numberArg
        # The node is void
        self.void = True
        # locate the node's neighbors
        self.neighbors = array(findNeighbors(numberArg, levelArg))
        # identify the node's father
        self.father = findFather(numberArg)
        # identify the node's childs
        self.childs = array(findChilds(numberArg))
        # center of the node
        n1, n2 = bitDeinterleaving(numberArg, levelArg)
        self.centerN1 = 2**(-levelArg)*(n1 + 0.5)
        self.centerN2 = 2**(-levelArg)*(n2 + 0.5)
        # find the scaled center of the node
        self.centerX = self.centerN1
        self.centerY = self.centerN2
        self.centerZ = complex(n1, n2)
        # Interaction list
        interactionL = array(findParentsChildren(numberArg, levelArg))
        self.interactionList = setdiff1d(interactionL, self.neighbors)
        # Create the blob arrays
        self.blobZ = array([], complex)
        self.blobVel = array([], complex)
        self.blobVelUpdate = array([], complex)
        self.blobCirculation = array([], float)
        self.blobNumber = array([])

    def __str__(self):
        s = 'num '+str(self.number)+'\n'
        s += 'center '+str(self.centerX)+','+str(self.centerY)+'\n'
        s += 'blob pos '+str(self.blobZ)+'\n'
        s += 'blob circ '+str(self.blobCirculation)+'\n'
        s += 'blob vel '+str(self.blobVel)
        return s
    
    def addBlob(self, z, vel, circulation, blobNum):
        '''Add information of a blob to the Quadtree'''
        self.blobZ = append(self.blobZ, z)
        self.blobVel = append(self.blobVel, vel)
        self.blobVelEval = append(self.blobVel, 0j)         # Saves dummy data in this variable
        self.blobCirculation = append(self.blobCirculation, circulation)
        self.blobNumber = append(self.blobNumber, blobNum)
        # The Node is not void anymore
        self.void = False
    
    def empty(self):
        ''' Empty the node parameters of the box'''
        self.blobZ = array([], complex)
        self.blobVel = array([], complex)
        self.blobVelUpdate = array([], complex)
        self.blobCirculation = array([], float)
        self.blobNumber = array([])

class Quadtree:
    '''
    Representation of a Quadtree, initialize the quadtree structure.
    '''
    def __init__(self, levelArg):
        '''
        Create the QuadtreeNodes nodes.
        '''
        # Contador de niveles comienza con nivel 0
        self.level = levelArg
        self.levels = [[],[]]
        # Create levels starting at level 2
        li = range(2, levelArg + 1)
        for levelNum in li:
            # Create the new nodes for the level levelNum
            # Add the new nodes to the level levelNum
            self.levels.append([QuadtreeNode(nodeNum, levelNum) for nodeNum in range(0, 4**levelNum)])

    def __str__(self):
        s = 'level '+str(self.level)+'\n'
        for box in self.levels[-1]:
            s += str(box)+'\n'
        return s

    def findNode(self, point, level):
        '''
        Find the number of the box at level \level\ where point \point\ belongs.
        Input:
            point - point consulted
            level - level of the quadtree being consulted
        '''
        nodeNumber = 0        
        # Scale the blobs to find Quadtree-Node
        xP = (point.real - self.minX) / self.diffX
        yP = (point.imag - self.minY) / self.diffY
        # apply interleaving and find the node at the corresponding level
        num1 = toBin(xP, level)
        num2 = toBin(yP, level)
        nodeNumber = bitInterleaving(num1, num2)
        return nodeNumber
            
    def fillQuadtree(self, arrZ, arrVel, arrCirculation):
        '''
        Fill the quadtree s nodes with blob information.
        Input:
            arrZ - complex array that contains the position of the blobs
            arrVel - complex array that contains the velocity of the blobs
            arrCirculation - float array that contains the circulation of each blob
        '''
        arrX = arrZ.real
        arrY = arrZ.imag
        self.N = len(arrX)
        # Finds min and max
        self.minX = min(arrX)
        self.maxX = max(arrX)
        self.diffX = self.maxX - self.minX
        self.minY = min(arrY)
        self.maxY = max(arrY)
        self.diffY = self.maxY - self.minY
        # Assign each blob to the associated Quadtree-Node
        for i in range(self.N):
            # Scale the blobs to find Quadtree-Node
            xP = (arrX[i] - self.minX) / self.diffX
            yP = (arrY[i] - self.minY) / self.diffY
            # apply interleaving and find the node at the corresponding level
            num1 = toBin(xP, self.level)
            num2 = toBin(yP, self.level)
            numInterleaved = bitInterleaving(num1, num2)
            #print 'x',xP,'y',yP,num1,num2
            #print 'blob',i,'box',numInterleaved
            # Add the blob to the Quadtree-node
            associatedNode = self.levels[self.level][numInterleaved]
            associatedNode.addBlob(arrZ[i], arrVel[i], arrCirculation[i], i)
        # Assign to each Node its not-scaled center coordinates
        for levelNum in range(2, self.level + 1): # From 2 to the lowest level
            liN = range(0, 4**levelNum)
            for nodeNum in liN:
                # recover the node
                quadtreeNode = self.levels[levelNum][nodeNum]
                # re-scale the node s center
                quadtreeNode.centerX = (self.diffX * quadtreeNode.centerX) + self.minX
                quadtreeNode.centerY = (self.diffY * quadtreeNode.centerY) + self.minY
                quadtreeNode.centerZ = complex(quadtreeNode.centerX, quadtreeNode.centerY)
        
    def emptyBoxes(self, boxesToEmpty):
        '''
        Void the boxes of the quadtree contained in the input array.
        Input:
            boxesToEmpty - array with the number id of the boxes to empty.
        '''
        maxlevel = self.levels[self.level]      # obtain the maxlevel of the QT
        for b in boxesToEmpty:
            node = maxlevel[b]
            node.empty()
## End File

