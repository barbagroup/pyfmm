"""
Filename: particleDistribution.py
Comment: Setup the initial position of particles
"""

from numpy import *
import random


def latticeDistribution(leftLimit, rightLimit, topLimit, bottomLimit, spacing):
    '''
    Obtains the two dimensional computational domain limits: 
        left, right, top, bottom
    and fill the computational domain with particles equispaced in
    the two dimensions.
    '''
    x , y = mgrid[leftLimit : rightLimit : spacing, topLimit : bottomLimit : spacing]
    x = x.flatten()
    y = y.flatten()
    
    z = []
    for i in range(len(x)):
        # Add the equispaced position to the particles' array
        z.extend([complex(x[i],y[i])])
        
    z = array(z)                # Convert the array to a numpy object
    return z


def triangleDistribution(leftLimit, rightLimit, topLimit, bottomLimit, spacing):
    '''
    Obtains the two dimensional computational domain limits: 
        left, right, top, bottom
    and fill the computational domain with particles equispaced in
    the two dimensions and then shift the even .
    '''
    halfSpacing = (spacing / 2)
    x1, y1 = mgrid[leftLimit: (rightLimit - halfSpacing) : spacing, topLimit : (bottomLimit - spacing) : 2 * spacing]
    
    x = x1.flatten()
    y = y1.flatten()
    
    size = len(x)
    z = []
    for i in range(size):
        # Add the equispaced position to the particles' array
        z.extend([complex(x[i],y[i])])
    
    
    x2, y2 = mgrid[(leftLimit + halfSpacing): rightLimit : spacing, (topLimit - spacing) : bottomLimit : 2 * spacing]
    x = x2.flatten()
    y = y2.flatten()
    
    size = len(x)
    for i in range(size):
        # Add the equispaced position to the particles' array
        z.extend([complex(x[i],y[i])])
        
    z = array(z)                # Convert the array to a numpy object
    return z



def quasiRandomDistribution(leftLimit, rightLimit, topLimit, bottomLimit, spacing, noise):
    '''
    Obtains the two dimensional computational domain limits: 
        left, right, top, bottom
    and fill the computational domain with particles quasi random equispaced in
    the two dimensions.
    '''
    x , y = mgrid[leftLimit : rightLimit : spacing, topLimit : bottomLimit : spacing]
    x = x.flatten()
    y = y.flatten()
    
    z = []
    halfSpacing = noise / 2
    
    for i in range(len(x)):
        # Generates uniform noise between '- half the spacing' to '+ half the spacing'
        #   in both dimensions
        xNoise = random.uniform(-halfSpacing, halfSpacing)
        yNoise = random.uniform(-halfSpacing, halfSpacing)
        # add the uniform noise to the particle position
        z.extend([complex(x[i] + xNoise, y[i] + yNoise)])
        
    z = array(z)                # Convert the array to a numpy object
    return z


# EOF