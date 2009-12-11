"""
Fast Multipole Method test.
Solving Vortex Blob Method for Lamb Oseen test case.
The number of particles of the problem depends on the
size of the computational domain

Usage:

-l <number>     Number of levels for the FMM

-p <number>     Truncation Number for the FMM

-n <number>     Size of the computational domain. Values 0 to 90

-d <number>     Particles distribution:
                0: Lattice distribution
                1: Triangular distribution
                2: Single particle distribution

-i <number>     Vorticity initialization
                0: Lamb Oseen
                1: Random vorticity
                2: Single particle

-h              Show help
"""

from numpy import *
import profile
import pstats
import os
import time
import getopt
import csv
import sys


## Import local modules  ----------------------------------------------
from pyFMM.FastMultipole.fastMultipoleMethod import FMMevalVelocity
from velocityBlobs import *
from simlogger import *
from particleDistribution import *


## Constants
# particles distribution
LATTICE_DIST = 0
TRIANGULAR_DIST = 1
QUASI_RAND_DIST = 2
# Vorticity initialization
LAMBOSEEN_INI = 0
RAMDOM_VORT_INI = 1
SINGLE_PART_INI = 2
# general constants
EPS = 10**(-20)     # epsilon machine
SIM_DOMAIN_SIZE = mgrid[1.:10.1:.1]

# function that computes the Lamb Oseen Vorticity
def lambOseen(gamma, nu, z, t):
    #centerLambOseen = 0.1 - 0.1j
    centerLambOseen = 0j
    r = abs(centerLambOseen - z)
    c0 = 4. * nu * t
    c1 = gamma / math.pi
    vort = (c1 / c0) * exp (- r**2 / c0)
    return vort

# function that computes the Lamb Oseen Velocity
def lambOseenVel(gamma, nu, z, t):
    centerLambOseen = 0j
    r = abs(z - centerLambOseen) + EPS
    nr2 = - r**2
    c0 = gamma / (2. * math.pi * r**2)
    c1 = nr2 / (4. * nu * t)
    vel = c0 * (1. - exp (c1))
    vel = (-z.imag + z.real * 1j) * vel
    return vel

def main():
    ## Default parameters
    ## FMM Parameters -----------------------------------------------
    level_param = 3
    p_param = 5
    h1_param = 0
    h2_param = 0
    # simulation parameters
    save_log = False                        # save sim data in the log file
    save_run = True                         # save sim information
    simulation_size = 0
    compare_analytical = 1                  # compare FMM vs analytical problem
    ## input program parameters -----------------------------------
    particle_distribution = LATTICE_DIST    # default distribution of the particles
    vorticity_distribution = LAMBOSEEN_INI  # default distribution of the initial vorticity
    gamma = 1.                              # gamma parameter of the lamb Oseen
    nu = 0.0005                             # viscosity value
    k = 2.                                  # blob parameter k
    sigma = 0.02                            # particle radius
    overlap = 0.8                           # overlap ratio
    tini = 4.                               # initial time for the simulation
    tend = 4.                               # end time of the simulation
    dt = 0.02                               # time step of the simulation
    steps = int((tend - tini) / dt)         # number of simulation steps
    s_grid = SIM_DOMAIN_SIZE[simulation_size]   # the side of the grid
    ## Variables calculations
    h = overlap * sigma                     # spacing of the particles
    sigma2 = sigma**2
    t = tini                                # time for the first step
    time_antiDiff = sigma2 / (2. * nu)      # Time "anti-diffusion" correction
    noise = 1. * h                         # noise value is X times the lattice step
    
    ## parse command line options------------------------------------
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'sghp:l:u:v:n:d:i:', ['help', 'saveGraph', 'level:'])
    except getopt.error, msg:
        print msg
        print "for help use --help"
        sys.exit(2)
    # process options
    for o, a in opts:
        # process help option
        if o in ('-h', '--help'):
            print __doc__
            sys.exit(0)
        # process truncation parameter
        if o in ('-p'):
            p_param = int(a)
        # process level parameter
        if o in ('-l', '--level'):
            level_param = int(a)
        # process FMM h1 shift parameter
        if o in ('-u'):
            h1_param = int(a)
        # process FMM h2 parameter
        if o in ('-v'):
            h2_param = int(a)
        # process simulation size parameter
        if o in ('-n'):
            simulation_size = int(a)
            if (simulation_size < len(SIM_DOMAIN_SIZE)):
                s_grid = SIM_DOMAIN_SIZE[simulation_size]
        # process particle distribution parameter
        if o in ('-d'):
            particle_distribution = int(a)
        # process vorticity distribution parameter
        if o in ('-i'):
            vorticity_distribution = int(a)
    
    ## Initialization -----------------------------------------------
    # create the grid that contains the coordinates x,y
    z = ''
    if particle_distribution == LATTICE_DIST:
        z = latticeDistribution(-s_grid/2, s_grid/2, -s_grid/2, s_grid/2, h)
        
    elif particle_distribution == TRIANGULAR_DIST:
        z = triangleDistribution(-s_grid/2, s_grid/2, -s_grid/2, s_grid/2, h)
        
    elif particle_distribution == QUASI_RAND_DIST:
        z = quasiRandomDistribution(-s_grid/2, s_grid/2, -s_grid/2, s_grid/2, h, noise)

    # Initialize the particles circulation
    w = lambOseen(gamma, nu, z, t - time_antiDiff)     # vorticity of the blob with time shifting fix
    circulation = w * h**2
    vel = zeros((len(z)), complex)
    
    print '\tNumber of blobs: ' + str(len(vel))
    print '\tTotal circulation in particles: ' + str(sum(circulation))
    
    ################## Experiment:
    
    timeTotal = time.clock()
    print 'FMM Started'
    # Calculate  velocity FMM
    fmmTime = time.clock()
    circulation, z, vel = FMMevalVelocity(level_param, p_param, circulation, z, vel, sigma2, k)
    fmmTime = time.clock() - fmmTime
    print 'FMM Finished'
    
    vel2 = []
    directTime = 0
    if (compare_analytical == 0):
        print 'Direct Started'
        # computation of the error against the DIRECT calculation
        directTime = time.clock()
        vel2 = evalVelocityArr(circulation, z, sigma2, k)
        directTime = time.clock() - directTime
        print 'Direct Finished'
    else:
        # computation of error agains the ANALYTIC solution
        vel2 = lambOseenVel(gamma, nu, z, t)
    
    error = vel - vel2
    errorRel = log10(EPS + abs(error) / (max(abs(vel2)) + EPS))
    errorL2 = linalg.norm(error) / linalg.norm(vel2)
    
    # computes the total time    
    timeTotal = time.clock() - timeTotal
    
    ############### End Experiment
    
    print 'Problem: LambOseen'
    print 'Vortex Blob Method:'
    print '\tNumber of blobs: ' + str(len(vel))
    print 'Fast Multipole Method:'
    print '\tNumber of levels: ' + str(level_param)
    print '\tTruncation number: ' + str(p_param)
    print '\tParticles per box: ' + str(len(vel)/(4**level_param))
    print '\tTotal circulation in particles: ' + str(sum(circulation))
    # Time Measuring
    print 'Time:'
    print '\t Direct: ' + str(directTime)
    print '\t FMM: ' + str(fmmTime)
    print '\t Total: ' + str(timeTotal)
    print 'Velocity:'
    print '\tMax vel: ' + str(max(abs(vel)))
    print '\tMin vel: ' + str(min(abs(vel)))
    print '\tMean vel: ' + str(mean(abs(vel)))
    print 'Relative Error:'
    print '\tLog Max error: ' + str(max(errorRel))
    print '\tLog Min error: ' + str(min(errorRel))
    print '\tError L2: ' + str(errorL2)
    
    maxError = -100
    maxErrorZ = 0j
    errorPosition = 0
    for i in range(len(errorRel)):
        if errorRel[i] > maxError:
            errorPosition = i
            maxError = errorRel[i]
            maxErrorZ = z[i]
    print 'Pos Max Error Rel: ' + str(maxErrorZ)
    print 'Direct Value:\t' + str(vel2[errorPosition])
    print 'FMM Value:\t' + str(vel[errorPosition])

    if save_run:
        sim_str = 'N' + str(len(vel)) + '_L' + str(level_param) + '_P' + str(p_param)
        print sim_str
        # Create a new Logger
        simulation_folder = 'Batch_' +'VORT' + str(vorticity_distribution) + '_DIST' + str(particle_distribution)
        logger = SimLogger('SimulationLogger' + '_N' + str(len(vel)), simulation_folder)
        # Save run info
        runInfo = logger.csvOutputLog('runData')
        runInfo.addElement(len(vel))                              # Number of Blobs
        runInfo.addElement(level_param)                           # FMM setup - Levels
        runInfo.addElement(p_param)                               # FMM setup - Truncation
        # Comparision against DIRECT calculation
        runInfo.addElement(max(errorRel))                         # save max velocity error
        runInfo.addElement(min(errorRel))                         # save min velocity error
        runInfo.addElement(errorL2)                               # save L-2 norm velocity error
        # Time results
        runInfo.addElement(fmmTime)                               # save time elapsed for FMM
        runInfo.addElement(directTime)                            # save time elapsed for DIRECT calculation
        runInfo.addElement(directTime/fmmTime)                    # acceleration ratio
        runInfo.flushRow()
        
        # Save error data
        logger.saveData(sim_str + '_simulation_Error', 
                POINTS_DIST_SCATTER, [z.real, z.imag, errorRel])               # Data from error
        

# Run Main function
if __name__ == "__main__":
    main()
