## Developing an open-source Python implementation of the Fast Multipole Algorithm (FMM) for scientific applications. ##

The FMM can be used in many scientific computing applications: the simulation of many stars, electrostatics, the calculation of atoms or molecules out of equilibrium, and particle methods for continuum problems. The advantages of the FMM can be huge when large numbers of particles are involved, as it reduces the complexity of calculations from O(N<sup>2</sup>) to O(N).

A more widespread adoption of the FMM algorithm has not occurred, mainly due to the complexity of the algorithm and the considerable programming effort, when compared with other algorithms like particle-mesh methods, or treecodes providing O(N log N) complexity.

We are developing an open source implementation of the FMM —with particular application to the calculation of a velocity field induced by N vortex particles.

**February 2008** -- at this time, we release a stable _alfa_ version of the Python code, and invite interested parties to email us if they would like to collaborate with us on further developments.

|We distribute this code under the MIT License, giving potential users the greatest freedom possible.  We do, however, request fellow scientists that if they use our codes in research, they kindly include us in the _acknowledgement_ of their papers.  We do not request gratuitous citations;  only cite our articles if you deem it warranted.|
|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|

### Related links ###

  * **[PetFMM](http://barbagroup.bu.edu/Barba_group/PetFMM.html)**: Open-source, parallel C++ implementation of the Fast Multipole Method. PetFMM can be obtained directly from its [repository](http://bitbucket.org/petfmm/petfmm-dev/) or from the [group's webpage](http://barbagroup.bu.edu/Barba_group/PetFMM.html).

### Publications ###

  * _"Characterization of the errors of the FMM in particle simulations"_  by Felipe A. Cruz and L. A. Barba.

> Preprint uploaded to [ArXiv](http://arxiv.org/abs/0809.1810) on 10 September 2008.

> Although the literature on the subject provides theoretical error bounds for the FMM approximation, there are not many reports of the measured errors in a suite of computational experiments.  We have performed such an experimental investigation, and summarized the results of about 1000 calculations using the FMM algorithm, to characterize the accuracy of the method in relation with the different parameters available to the user.  In addition to the more standard diagnostic of the maximum error, we supply illustrations of the spatial distribution of the errors, which offers visual evidence of all the contributing factors to the overall approximation accuracy: multipole expansion, local expansion, hierarchical spatial decomposition (interaction lists, local domain, far domain). This presentation is a contribution to any researcher wishing to incorporate the FMM acceleration to their application code, as it aids in understanding where accuracy is gained or compromised.

  * _"Characterization of the accuracy of the fast multipole method in particle simulations"_

> This is the published paper, a revised version of the ArXiV preprint linked above.

> Published online: May 5, 2009.  [DOI](http://dx.doi.org/10.1002/nme.2611)