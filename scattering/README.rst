scattering - single scattering properties of particles
#########


This is the python module scattering.
This module is part of the larger pyPamtra2 software suite for the simulation of consistent radiative transfer between active and passive observations mainly at microwave frequencies.
The module computes the single scattering properties of hydrometeor models using diverse electromagnetic scattering simulators.
The functionalities of this module can be used outside of pyPamtra2 in any other application.

The general idea of the module is that it is possible to construct an abstraction layer between radiative transfer codes and single
scattering models. Regardless of the single scattering solution employed the resulting products can be geralized in terms of cross sections
and angle resolved scattering matrices. This module implements this level of abstraction and attempts to increase efficiency of the code
by selectively avoiding the calculation of certain quantities that are unable to be represented by certain scattering approximations

include
=======

* Rayleigh approximation for spheres
* Rayleigh approximation for spheroids
* Lorenz-Mie scattering solution for spheres
* Self-Similar Rayleigh Gans
* T-Matrix for spheroids
* Scattering DataBases (DBs) - Liu, Hong, Leinonen


