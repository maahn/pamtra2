#! /usr/bin/env bash 
set -e

cd dfftpack 
make  
cd ../

cd pyPamtraRadarSimulator
python setup.py install
cd ../

cd pyPamtraRadarMoments
python setup.py install
cd ../

cd pyPamtra2
python setup.py install
cd ../

cd refractive
python setup.py install
cd ../

cd scattering
python setup.py install
cd ../
