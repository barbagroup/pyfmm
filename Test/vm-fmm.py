#!/usr/bin/env python
'''
Velocity computations with Vortex method
using the Fast Multipole Method 
'''

from numpy import *

## Import local modules
from pyFMM.DataStructure.quadtree import *
from pyFMM.FastMultipole.fastMultipoleMethod import FMMevalVelocity

import time
import os, sys
sys.path.insert(0, os.path.join(os.environ['PETSC_DIR'], 'config', 'BuildSystem'))
import script
import RDict

# general constants
EPS = 10**(-20)     # epsilon machine

class Options(script.Script):
    def __init__(self):
        script.Script.__init__(self, clArgs = ['-log=verify.log']+sys.argv[1:], argDB = RDict.RDict())
        self.setup()
        return

    def setupHelp(self, help):
        '''This method should be overidden to provide help for arguments'''
        import nargs

        script.Script.setupHelp(self, help)
        help.addArgument('FMM', '-debug',    nargs.ArgBool(None, False, 'Debug FMM', isTemporary = 1))
        help.addArgument('FMM', '-lower',    nargs.Arg(None, [-0.5, -0.5],  'Lower-left corner of the tree', isTemporary = 1))
        help.addArgument('FMM', '-upper',    nargs.Arg(None, [ 0.5,  0.5],  'Upper-right corner of the tree', isTemporary = 1))
        help.addArgument('FMM', '-maxLevel', nargs.ArgInt(None, 2,  'Last tree level (numLevels-1)', 0, isTemporary = 1))
        help.addArgument('FMM', '-terms',    nargs.ArgInt(None, 5,  'Number of terms in multipole expansion', 1, isTemporary = 1))
        help.addArgument('FMM', '-sigma',    nargs.ArgReal(None, 0.02,  'Particle radius', 0.0, isTemporary = 1))
        help.addArgument('FMM', '-k',        nargs.ArgInt(None, 2,  'Blob parameter???', 0, isTemporary = 1))
        help.addArgument('FMM', '-overlap',  nargs.ArgReal(None, 0.8,  'Overlap ratio', 0.0, isTemporary = 1))
        help.addArgument('FMM', '-gamma',    nargs.ArgReal(None, 1.0,  'Gamma parameter of the Lamb-Oseen model', 0.0, isTemporary = 1))
        help.addArgument('FMM', '-nu',       nargs.ArgReal(None, 0.0005,  'Viscosity parameter of the Lamb-Oseen model', 0.0, isTemporary = 1))
        help.addArgument('FMM', '-t',        nargs.ArgReal(None, 4.0,  'Time at first step', isTemporary = 1))
        help.addArgument('FMM', '-true_solution_type', nargs.Arg(None, 'analytic',  'True solution type: direct, analytic', isTemporary = 1))
        help.addArgument('FMM', '-basename', nargs.Arg(None, '',  'Base filename for output', isTemporary = 1))
        return help

    def getDebug(self):
        if not hasattr(self, '_debug'):
            self._debug = self.argDB['debug']
        return self._debug
    def setDebug(self, debug):
        self._debug = debug
        return
    debug = property(getDebug, setDebug, doc = 'Debug FMM')
    def getLower(self):
        if not hasattr(self, '_lower'):
            self._lower = map(float, self.argDB['lower'])
        return self._lower
    def setLower(self, lower):
        self._lower = lower
        return
    lower = property(getLower, setLower, doc = 'Lower-left corner of the tree')
    def getUpper(self):
        if not hasattr(self, '_upper'):
            self._upper = map(float, self.argDB['upper'])
        return self._upper
    def setUpper(self, upper):
        self._upper = upper
        return
    upper = property(getUpper, setUpper, doc = 'Upper-right corner of the tree')
    def getMaxLevel(self):
        if not hasattr(self, '_maxLevel'):
            self._maxLevel = self.argDB['maxLevel']
        return self._maxLevel
    def setMaxLevel(self, maxLevel):
        self._maxLevel = maxLevel
        return
    maxLevel = property(getMaxLevel, setMaxLevel, doc = 'Last tree level (numLevels-1)')
    def getTerms(self):
        if not hasattr(self, 'terms'):
            self._terms = self.argDB['terms']
        return self._terms
    def setTerms(self, terms):
        self._terms = terms
        return
    terms = property(getTerms, setTerms, doc = '')
    def getSigma(self):
        if not hasattr(self, '_sigma'):
            self._sigma = self.argDB['sigma']
        return self._sigma
    def setSigma(self, sigma):
        self._sigma = sigma
        return
    sigma = property(getSigma, setSigma, doc = '')
    def getK(self):
        if not hasattr(self, '_k'):
            self._k = self.argDB['k']
        return self._k
    def setK(self, k):
        self._k = k
        return
    k = property(getK, setK, doc = '')
    def getOverlap(self):
        if not hasattr(self, '_overlap'):
            self._overlap = self.argDB['overlap']
        return self._overlap
    def setOverlap(self, overlap):
        self._overlap = overlap
        return
    overlap = property(getOverlap, setOverlap, doc = '')
    def getGamma(self):
        if not hasattr(self, '_gamma'):
            self._gamma = self.argDB['gamma']
        return self._gamma
    def setGamma(self, gamma):
        self._gamma = gamma
        return
    gamma = property(getGamma, setGamma, doc = '')
    def getNu(self):
        if not hasattr(self, '_nu'):
            self._nu = self.argDB['nu']
        return self._nu
    def setNu(self, nu):
        self._nu = nu
        return
    nu = property(getNu, setNu, doc = '')
    def getT(self):
        if not hasattr(self, '_t'):
            self._t = self.argDB['t']
        return self._t
    def setT(self, t):
        self._t = t
        return
    t = property(getT, setT, doc = '')
    def getTrueSoln(self):
        if not hasattr(self, '_trueSoln'):
            self._trueSoln = self.argDB['true_solution_type']
        return self._trueSoln
    def setTrueSoln(self, trueSoln):
        self._trueSoln = trueSoln
        return
    trueSoln = property(getTrueSoln, setTrueSoln, doc = '')
    def getBasename(self):
        if not hasattr(self, '_basename'):
            self._basename = self.argDB['basename']
        return self._basename
    def setBasename(self, basename):
        self._basename = basename
        return
    basename = property(getBasename, setBasename, doc = '')

    def __str__(self):
        s = ['FMM Options:',
             '  lower:     '+str(self.lower),
             '  upper:     '+str(self.upper),
             '  numLevels: '+str(self.maxLevel+1),
             '  p:         '+str(self.terms),
             '  sigma:     '+str(self.sigma),
             '  k:         '+str(self.k),
             '  overlap:   '+str(self.overlap),
             '  gamma:     '+str(self.gamma),
             '  nu:        '+str(self.nu),
             '  t:         '+str(self.t),
             '  trueSoln:  '+str(self.trueSoln),
             '  basename:  '+str(self.basename)]
        return '\n'.join(s)

def evalVelocityArr(circulation, Z, sigma2, k):
    '''
    Direct evaluation of the velocity from the vortex blobs
    by using the Biot-Savart Law.
    '''
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


def lambOseen(gamma, nu, z, t):
    ''' Computes the Lamb-Oseen vortex vorticity centered at origin.
    gamma   - float Lamb-Oseen vortex parameter.
    nu      - float Lamb-Oseen vortex parameter.
    z       - complex evaluation point.
    t       - float evaluation time.
    '''
    centerLambOseen = 0j
    r = abs(centerLambOseen - z)
    c0 = 4. * nu * t
    c1 = gamma / math.pi
    vort = (c1 / c0) * exp (- r**2 / c0)
    return vort

# function that computes the Lamb Oseen Velocity
def lambOseenVel(gamma, nu, z, t):
    centerLambOseen = 0j
    r = abs(z - centerLambOseen)
    nr2 = - r**2
    c0 = gamma / (2. * math.pi * (r**2 + EPS))
    c1 = nr2 / (4. * nu * t)
    vel = c0 * (1. - exp (c1))
    vel = (-z.imag + z.real * 1j) * vel
    return vel
    
def latticeDistribution(leftLimit, rightLimit, bottomLimit, topLimit, spacing):
    '''
    Obtains the two dimensional computational domain limits: 
        left, right, top, bottom
    and fill the computational domain with particles equispaced in
    the two dimensions.
    '''
    x , y = mgrid[leftLimit+spacing/2.0 : rightLimit : spacing, bottomLimit+spacing/2.0 : topLimit : spacing]
    x = x.flatten()
    y = y.flatten()
    
    z = []
    for i in range(len(x)):
        # Add the equispaced position to the particles' array
        z.extend([complex(x[i],y[i])])
        
    z = array(z)                # Convert the array to a numpy object
    return z

def verify(options, z, tree, vel, velDirect):
    '''Output a verification file'''
    f = open(options.basename+'fmm.'+str(len(vel))+'.verify', 'w')
    # Number of levels
    f.write(str(options.maxLevel+1)+'\n')    
    # Number of terms
    f.write(str(options.terms)+'\n')    
    # Number of blobs
    f.write(str(len(vel))+'\n')
    # Blob boxes
    x = z.real
    y = z.imag
    for i in range(tree.N):
        # Scale the blobs to find Quadtree-Node
        xP = (x[i] - tree.minX) / tree.diffX
        yP = (y[i] - tree.minY) / tree.diffY
        # apply interleaving and find the node at the corresponding level
        num1 = toBin(xP, tree.level)
        num2 = toBin(yP, tree.level)
        numInterleaved = bitInterleaving(num1, num2)
        f.write(str(i)+' '+str(numInterleaved)+'\n')
    # Quadtree
    f.write(str(tree.minX)+'\n'+str(tree.maxX)+'\n'+str(tree.diffX)+'\n')
    f.write(str(tree.minY)+'\n'+str(tree.maxY)+'\n'+str(tree.diffY)+'\n')
    for l in range(options.maxLevel+1):
        f.write(str(l)+'\n')
        for box in tree.levels[l]:
            f.write(str(box.number)+'\n')
            f.write(str(box.centerX)+' '+str(box.centerY)+'\n')
            f.write(str(len(box.blobZ))+'\n')
            f.write(' '.join([str(c) for c in box.childs])+'\n')
            f.write(' '.join([str(n) for n in box.neighbors])+'\n')
            f.write(' '.join([str(i) for i in box.interactionList])+'\n')
            f.write(' '.join(['('+str(c.real)+','+str(c.imag)+')' for c in box.C])+'\n')
            f.write(' '.join(['('+str(d.real)+','+str(d.imag)+')' for d in box.D])+'\n')
            f.write(' '.join(['('+str(e.real)+','+str(e.imag)+')' for e in box.E])+'\n')
    # Direct solution
    for vD in velDirect:
        f.write('('+str(vD.real)+','+str(vD.imag)+')\n')
    # FMM solution
    for v in vel:
        f.write('('+str(v.real) +','+str(v.imag) +')\n')
    # log(||e||_2 / ||v^*||_2)
    error   = vel - velDirect
    errorL2 = log10(linalg.norm(error) / linalg.norm(velDirect))
    f.write(str(errorL2)+'\n')
    f.close()
    return
    
def main():
    options = Options()
    # Calculated variables
    h = options.overlap * options.sigma       # spacing of the particles
    sigma2 = options.sigma**2
    time_antiDiff = sigma2 / (2. * options.nu)      # Time "anti-diffusion" correction
    
    # Initialize particles positions
    z = latticeDistribution(options.lower[0], options.upper[0], options.lower[1], options.upper[1], h)
    ## Initialize the particles circulation
    # vorticity of the blob with time shifting fix (improves initialization)
    ### THIS SHOULD BE INITIALIZED BY ANALYTIC IF YOU ARE USING IT
    w = lambOseen(options.gamma, options.nu, z, options.t - time_antiDiff)
    circulation = w * h**2
    vel = zeros((len(z)), complex)
    
    # Calculate initial velocity
    fmmTime = time.clock()
    circulation, z, tree, vel = FMMevalVelocity(options.maxLevel, options.terms, circulation, z, vel, sigma2, options.k)
    fmmTime = time.clock() - fmmTime
    
    # computation of the error against the true solution
    directTime = time.clock()
    if options.trueSoln == 'analytic':
      true_vel = lambOseenVel(options.gamma, options.nu, z, options.t)
    elif options.trueSoln == 'direct':
      true_vel = evalVelocityArr(circulation, z, sigma2, options.k)
    directTime = time.clock() - directTime
    
    # FMM approximation relative error
    ##print 'vel',vel
    ##print 'vel direct',true_vel
    error = vel - true_vel
    error_rel = log10(EPS + abs((error) / (max(abs(true_vel)) + EPS)))
    errorL2 = log10(linalg.norm(error) / linalg.norm(true_vel))

    print options
    print 'Number of blobs:       ',len(vel)
    print 'FMM Relative Error L2: ',errorL2
    print 'FMM time:              ',fmmTime
    print 'Direct time:           ',directTime
    verify(options, z, tree, vel, true_vel)
    return

# Run Main function
if __name__ == "__main__":
    main()

## EOF
