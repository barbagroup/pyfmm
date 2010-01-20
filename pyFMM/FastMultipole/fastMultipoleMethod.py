"""
Implementacion of Upward and Downward pass. With combinatorics cache used.
"""

from numpy import *
from pyFMM.Kernels.translations import *
from pyFMM.DataStructure.quadtree import *
import time


def upwardStep1(quadtree, p, k):
    '''
    Upward pass Step 1:
    Generation of the multipole expansion coefficients for
    each box at the lowest level L
    Input:
        quadtree - quadtree structure that contains all the blobs ordered
        p - truncation number
        k - parameter of the blobs
    '''
    # Lowest level
    L = quadtree.level
    LowestLevel = quadtree.levels[L]
    # Loop over all the nodes of the lowest level
    for nodeB in LowestLevel:
        # get the center of the node
        Zc = nodeB.centerZ
        # C = 0
        C = zeros((p), complex)
        # Constant of the vortex blobs (simplified to points... explaination on paper)
        constant = -1j / (k * math.pi)
        # Loop over all the blobs contained on nodeB
        for i in range(len(nodeB.blobZ)):
            # blob position
            Zi = nodeB.blobZ[i]
            # get the multipole expansion coef(MEC)
            mec = MEC(p, Zi - Zc)
            # add the contribution of this blob to the total
            C = C + constant * dot(nodeB.blobCirculation[i], mec)
        # Save the Multipole Expansion coefficients of the nodeB
        nodeB.C = C.copy()
    return
    
    
    
def upwardStep2(quadtree, p):
    '''
    Upward pass step 2:
    Generation of the multipole expansions for every box at every level
    starting at level L-1 to level 2 calculating the re-expansions coef
    from the children s expansions in a recursive way.
    Input:
        quadtree - quadtree structure that contains all the blobs ordered
        p - truncation number
    '''
    # Lowest level
    L = quadtree.level
    # Loop levels starting at L-1 to 2
    levelsLoop = L + 1 - array(range(2,L))
    for l in levelsLoop:
        # loop over the nodes of the level
        for nodeB in quadtree.levels[l]:
            # get the center of the node
            Zc = nodeB.centerZ
            # Cnl = 0
            C = zeros((p), complex)
            # Loop over nodeB s children
            for childB in nodeB.childs:
                # get the children
                child = quadtree.levels[l+1][childB]
                # get the children center
                Zc_child = child.centerZ
                # calculate the multipole to multipole (M2M) translation matrix
                m2m = M2M(p, Zc - Zc_child)
                # add the contribution of this child to the total
                C = C + dot(m2m, child.C)
            # Save the Multipole Expansion coefficients of the nodeB
            nodeB.C = C
    
    return


def downwardStep1(quadtree, p):
    '''
    Downward Pass Step 1:
    Conversion of a multpole expansion into a local expansion
    Input:
        quadtree - quadtree structure that contains all the blobs ordered
        p - truncation number
    '''
    # Lowest level
    L = quadtree.level
    # Loop levels starting at 2 to L
    levelsLoop = array(range(2, L + 1))
    for l in levelsLoop:
        # loop over the nodes of the level
        for nodeB in quadtree.levels[l]:
            # get the center of the node
            Zc = nodeB.centerZ
            # Dnl = 0
            D = zeros((p), complex)
            # Loop over nodeB s interaction list
            for interactor in nodeB.interactionList:
                # get the interactor
                inter = quadtree.levels[l][interactor]
                # get the interactor s center
                Zc_inter = inter.centerZ
                # calculate the multipole to local (M2L) translation matrix
                m2l = M2L(p, Zc - Zc_inter)
                # add the contribution of this interactor to the total
                D = D + dot(m2l, inter.C)
            # Save the Multipole Expansion coefficients of the nodeB
            nodeB.D = D
    return


def downwardStep2(quadtree, p):
    '''
    Downward Pass Step 2:
    Translation of the local expansions.
    Input:
        quadtree - quadtree structure that contains all the blobs ordered
        p - truncation number
    '''
    # Lowest level
    L = quadtree.level
    # Loop levels starting at 2 to L
    levelsLoop = array(range(2, L + 1))
    for l in levelsLoop:
        # loop over the nodes of the level
        for nodeB in quadtree.levels[l]:
            if l == 2:
                # nodes in level 2 doesn't have parents... skip step
                nodeB.E = nodeB.D
            else:
                # get the center of the current node
                Zc = nodeB.centerZ
                # get the parent node
                parent = quadtree.levels[l-1][nodeB.father]
                # get the center of the parent node
                Zc_parent = parent.centerZ
                # calculate the local to local (L2L) translation matrix
                l2l = L2L(p, Zc - Zc_parent)
                # add the parent s contribution and nodeB s contributions to the total
                nodeB.E = nodeB.D + dot(l2l, parent.E)
    return


def evalVelocity2(circulation, Z, sigma2, k):
    # local variables
    size = len(Z)
    velocity = zeros((size), complex)
    # computation of common factors
    c1 = -1/(k * sigma2)
    c2 = 1j/(k * math.pi)
    for i in range(size):
        r2 = abs(Z - Z[i])**2
        r2[i] = 2**(-40)
        zz = (Z - Z[i]).conj()
        zz[i] = 1       # fix the division by zero problem
        # Calculates the velocity induced by the blob i of the center box
        blobInfluence = circulation[i] * c2 * (1 - exp(r2  * c1)) / zz
        blobInfluence[i] = 0    # blob self influence = 0
        # Calcules the velocity
        velocity = velocity + blobInfluence
    return velocity



def evalVelocity(circulationNB, ZNB, circulationCB, ZCB, sigma2, k):
    '''
    Calculates the influence of the particles in the neighborhood and center/self box.
    Input:
        circulationNB - array with the circulation of the blobs in the neighborhood (NB)
        ZNB - position of the blob in the neighborhood (NB) in complex representation
        circulationCB - array with the circulation of the blobs in the central box (CB)
        ZCB - position of the blob in the central box (CB) in complex representation
        sigma2 - sigma square parameter of the gaussian blob
        k - parameter k of the Gaussian blob
    Output:
        velocity - an array containing the velociy at each evaluation point ZCB
    '''
    # local variables
    sizeNB = len(ZNB)
    sizeCB = len(ZCB)
    velocity = zeros((sizeCB), complex)
    # computation of common factors
    c1 = -1/(k * sigma2)
    c2 = 1j/(k * math.pi)
    ###########################################
    # computation over the central box
    ###########################################
    for i in range(sizeCB):
        r2 = abs(ZCB - ZCB[i])**2
        r2[i] = 2**(-40)
        zz = (ZCB - ZCB[i]).conj()
        zz[i] = 1       # fix the division by zero problem
        # Calculates the velocity induced by the blob i of the center box
        blobInfluence = circulationCB[i] * c2 * (1 - exp(r2  * c1)) / zz
        blobInfluence[i] = 0    # blob self influence = 0
        # Calcules the velocity
        velocity = velocity + blobInfluence
    ##########################################
    # computation over the local neighborhood
    ##########################################
    for i in range(sizeNB):
        r2 = abs(ZCB - ZNB[i])**2
        zz = (ZCB - ZNB[i]).conj()
        # Calculates the velocity induced by the blob i of the neighborhood
        blobInfluence = circulationNB[i] * c2 * (1 - exp(r2  * c1)) / zz
        # add it to the induced velocity
        velocity = velocity + blobInfluence
    return velocity


def finalStep(quadtree, p, sigma2, k):
    '''
    Final Step Evaluation:
    Evaluation at all the points of the quadtree.
    Input:
        quadtree - quadtree structure that contains all the blobs ordered
        p - truncation number
        sigma2 - sigma square parameter of the gaussian blob
        k - parameter k of the Gaussian blob    
    '''
    # Lowest level
    L = quadtree.level
    # loop over the nodes of the level L
    for nodeB in quadtree.levels[L]:
        # Arrays that saves the blobs data for the center box (CB)
        ZCB = []
        CirculationCB = []
        # Arrays that saves the blobs data for the Neighborhood (NB)
        ZNB = []
        CirculationNB = []
        # Get the local blobs and circulation
        ZCB.extend(nodeB.blobZ)
        CirculationCB.extend(nodeB.blobCirculation)
        ###############################################
        # DIRECT EVALUATION PART
        ###############################################
        # Get the neighborhood nodes
        nodesNB = nodeB.neighbors
        # Get the neighborhood blobs and cirulation
        for nodeNum in nodesNB:
            # get the actual node
            node = quadtree.levels[L][nodeNum]
            ZNB.extend(node.blobZ)
            CirculationNB.extend(node.blobCirculation)
        # Calculate the direct interactions
        directVelocity = evalVelocity(CirculationNB, ZNB, CirculationCB, ZCB, sigma2, k)
        ###############################################
        # FMM EVALUATION PART
        ###############################################
        # Get the center of the box
        Zc = nodeB.centerZ
        fmmVelocity = zeros((len(ZCB)), complex)
        # loop over all the particles contained on the center box
        for i in range(len(ZCB)):
            local = Local(p, ZCB[i] - Zc)
            # Calculates the \conjugate velocity\ with the FMM and then in conjugate it
            # to obtain the \velocity\
            fmmVelocity[i] = (dot(nodeB.E, local)).conjugate()
        ###############################################
        # ADD EVERYTHING AND SAVE IT
        ###############################################
        nodeB.blobVelUpdate = directVelocity + fmmVelocity
    return


def FMMevalVelocity(level_param, p_param, circulation, z, vel, sigma2, k):
    '''
    Evaluates the velocity using the Vortex Blob Method accelerated with
    the Fast Multipole Method.
    Input:
        level_param - quadtree level. Parameter of the FMM
        p_param - truncation number. Parameter of the FMM.
        circulation - array with the circulation information of the blobs
        z - array with the positional information of the blobs
        vel - array with the velocity information of the blobs
        sigma2 - sigma square parameter of the gaussian blob
        k - parameter k of the Gaussian blob
    '''
    initCombCache(2*p_param)
    tree = Quadtree(level_param)
    # Filling the quadtree
    tree.fillQuadtree(z, vel, circulation)
    ##print tree
    ######## TEST ########
    #tree.empty([0])
    ######## END TEST ########
    # Start the FMM procedure
    upwardStep1(tree, p_param, k)
    upwardStep2(tree, p_param)
    downwardStep1(tree, p_param)
    downwardStep2(tree, p_param)
    finalStep(tree, p_param, sigma2, k)
    # Reconstruct the arrays
    out_number = []
    out_z = []
    out_vel = []
    out_circulation = []
    for nodeNum in range(4**level_param):
        nodeB = tree.levels[level_param][nodeNum]
        out_number.extend(nodeB.blobNumber)
        out_z.extend(nodeB.blobZ)
        out_vel.extend(nodeB.blobVelUpdate)
        out_circulation.extend(nodeB.blobCirculation)
    # Ordering the final output
    out_matrix = array([out_number,out_circulation, out_z, out_vel])
    out_matrix = out_matrix.transpose()
    out_matrix = out_matrix[out_matrix[:,0].argsort(),]
    out_circulation = out_matrix[:,1]
    out_z = out_matrix[:,2]
    out_vel = out_matrix[:,3]
    return out_circulation.real, out_z, tree, out_vel

    
