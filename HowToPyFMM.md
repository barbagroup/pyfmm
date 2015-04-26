This is  a short guide on how to install and use python FMM.

# Introduction #

Before downloading pyFMM check that your system fulfill all the requirements. This brief guide will introduce you on how to install, and use pyFMM.

## Pre-requisites ##

  * Python 2.5
  * Numpy
  * Scipy (recommended)


## Installation ##

Installation procedure:

### Linux / Mac OS X ###

To download the latest source from the svn (you need to have installed svn), from a terminal run:
```
svn checkout http://pyfmm.googlecode.com/svn/trunk/ pyfmm-read-only
```
As an alternative, you can download the source of the current version v1.X.X as a tar.gz file. Once downloaded, you can decompress the source from a terminal:
```
tar -xzf pyFMMv1.X.X.tar.gz
```

After downloading (and extracting) the latest version of the code, you need to add the folder of the package to the PYTHONPATH system variable, e.g.:
Let suppose that the final folder is /Users/user/pyFMMv1.X.X
```
export PYTHONPATH=${PYTHONPATH}:/Users/user/pyFMMv1.X.X
```
After this, the python interpreter should be able to detect the pyFMM package. So you are ready to test the installation. To do this, you can run our test program that is located in the Test folder of the package:
```
python pyFMMv1.X.X/Test/LambOseenTest.py
```
After a few seconds you should see the following output:
```
        Number of blobs: 3969
        Total circulation in particles: 1.0
FMM Started
FMM Finished
Problem: LambOseen
Vortex Blob Method:
        Number of blobs: 3969
Fast Multipole Method:
        Number of levels: 3
        Truncation number: 5
        Particles per box: 62
        Total circulation in particles: 1.0
Time:
         Direct: 0
         FMM: 5.4
         Total: 5.4
Velocity:
        Max vel: 1.13565329313
        Min vel: 0.112533052315
        Mean vel: 0.478659772956
Relative Error:
        Log Max error: -2.57063153267
        Log Min error: -5.75810088547
        Error L2: 0.000630939850659
Pos Max Error Rel: (-0.004+0.252j)
Direct Value:   (-0.631183237873-0.0100187815535j)
FMM Value:      (-0.634161869129-0.0106838376279j)
N3969_L3_P5
```