refractive - collection of methods to compute refractive index of hydrometeors
#########


This is the python module refractive.
This module is part of the larger pyPamtra2 software suite for the simulation of consistent radiative transfer between active and passive observations mainly at microwave frequencies.
The module computes the complex refractive index of ice and water between 1 and 1000 GHz.
The functionalities of this module can be used outside of pyPamtra2 in any other application.

Ice refractive index methids include:
====================================

* Matzler 2006 - Thermal Microwave Radiation: application to remote sensing
* Warren and Brandt 2008 - 'Optical constants of ice from the ultraviolet to the microwave: A revised compilation.' J. Geophys. Res., 113, D14220, doi:10.1029/2007JD009744. which updates and corrects Warren, S. G. (1984), 'Optical constants of ice from the ultraviolet to the microwave', Appl. Opt., 23, 1206–1225.
* Iwabuchi and Yang 2011 - 'Temperature dependence of ice optical constants: Implications for simulating the single-scattering properties of cold ice clouds' J. Quant. Spec. Rad. Tran. 112, 2520-2525

Water refractive index methods include:
=======================================

* Ellison (year) - biblio

Effective Medium Approximations (EMAs) include:
==============================================

* Maxwell-Garnett
* Bruggeman
* Sihvola

The module also includes utilities to easily compute effective dielectric properties of snow according to its effective density by assuming a homogeneous mixture of ice and air

The utilities.py submodule contains function definitions for dielectric properties manipulation

