"""
Calculates the velocity for the blobs

Input: 
    circulation - array with the circulation of the blobs
    Z - position of the blob in complex representation
    sigma2 - sigma square parameter of the gaussian blob
    k - parameter k of the Gaussian blob
Output:
    velocity - an array containing the velociy at each evaluation point
"""
from numpy import *

def evalVelocityArr(circulation, Z, sigma2, k):
    # local variables
    size = len(Z)
    velocity = zeros((size), complex)
    
    # common factors
    c1 = -1/(k * sigma2)
    c2 = 1j/(k * math.pi)

    for i in range(size):
        r2 = abs(Z - Z[i])**2
        r2[i] = 2**(-40)
        zz = (Z - Z[i]).conj()
        zz[i] = 1       # fix the division by zero problem
        # Calculates the velocity induced by the blob i
        blobInfluence = circulation[i] * c2 * (1 - exp(r2  * c1)) / zz
        blobInfluence[i] = 0    # blob self influence = 0
        # Calcules the velocity
        velocity = velocity + blobInfluence
    return velocity
    