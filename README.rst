Pamtra2
#######

EARLY ALPHA VERSION! NOT FEATURE COMPLETE!

Included Modules
================

* pyPamtra2 - The core package.
* pyPamtraRadarSimulator - Simulate a radar Doppler spectrum
* pyPamtraRadarMoments - Estimate the moments of the radar Doppler spectrum
* refractive - general purpose refractive index of ice, water and Effective Medium Approximations
* scattering - general purpose electromagnetic single-scattering module. Includes several scattering methods

Installation
============

run install.sh or use the setup.py scripts in the project folders

Single package installation
---------------------------
Change directory to the interested package and run either python or python3
```
python setup.py install [flag]
```

the flag --user will install the package under user's home and it is useful if you do not have root priviliges

Whole pyPamtra suite installation
---------------------------------
The install.sh script is provided for an easy installation of the whole packages distributed with pyPamtra
By default it uses python and installs system wide. Just run
```
./install.sh
```
You might want to pass some options to the script enabling python3 and/or user installation. To do so, run
```
./install.sh -p python3 -f --user
```
