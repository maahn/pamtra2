#! /usr/bin/env bash 
set -e


cd pyPamtraRadarSimulator
python setup.py install
cd ../

cd pyPamtraRadarMoments
python setup.py install
cd ../

cd pyPamtra2
python setup.py install
cd ../

cd dfftpack 
make  
cd ../

cd refractive
make
cd ../
