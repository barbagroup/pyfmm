"""
Implementation operators:
- Implementation of the Local expansion
- Implementation of the multipole expansions coefficients
- Implementation of the translation local to local
- Implementation of the translation Multipole to local
- Implementation of the translation Multipole to multipole
"""
from numpy import *

####################### support functions ##########################
# factorial function
fac = lambda n: [1, 0] [n > 0] or fac(n - 1) * n

# combinatorial function
def comb2(n, k):
    return (fac(n) / (fac(n - k) * fac(k)))

combCache = []

def initCombCache(combCacheSize):
    global combCache
    iSize = 2 * combCacheSize
    jSize = combCacheSize
    combCache = zeros((iSize, jSize), float)
    for i in range(iSize):
        for j in range(jSize):
            combCache[i][j] = comb2(i,j)
    return
    
def comb(n,k):
    return combCache[n][k]

####################################################################


####################################################################
#################### MULTIPOLE EXPANSION COEF ######################
####################################################################

def MEC(p, t):
    '''
    Computation of p multipole expansion coefficients
    Input:
        p - number of coefficients to compute
        t - translation distance
    ''' 
    mec = zeros((p), complex)
    for m in range(p):
        mec[m] = t**m
    return mec

####################################################################
###################### MULTIPOLE TO MULTIPOLE ######################
####################################################################

def M2M(p, t):
    '''
    Translation multipole to multipole
    Input: 
        p - number of p data
        t - translation distance
    '''
    translation = zeros((p,p), complex)
    for n in range(p):
        for m in range(n + 1):
           #translation[n][m] = (-t)**(n - m) * (comb(n, m))
           translation[n][m] = (-t)**(n - m) * (combCache[n][m])
    return translation
    

####################################################################
###################### MULTIPOLE TO LOCAL ##########################
####################################################################

def M2L(p, t):
    '''
    Translation multipole to local
    Input: 
        p - number of p data
        t - translation distance
    '''
    translation = zeros((p,p), complex)
    for n in range(p):
        for m in range(p):
            translation[n][m] = ((-1)**n) * (combCache[n + m][m]) / (t)**(m + n + 1)
    return translation

####################################################################
###################### LOCAL TO LOCAL ##############################
####################################################################

def L2L(p, t):
    '''
    Translation local to local
    Input: 
        p - number of p data
        t - translation distance
    '''
    translation = zeros((p,p), complex)
    for n in range(p):
        for m in range(n, p):
            translation[n][m] = (t**(m - n)) * (combCache[m][n])
    return translation

####################################################################
###################### LOCAL EXPANSION #############################
####################################################################

def Local(p, z):
    '''
    Computation of the Local expansion
    Input: 
        p - number of p data
        t - translation distance
    '''
    local = ones((p), complex)
    for i in range(1,p):
        local[i] = z**i
    return local


