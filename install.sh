#! /usr/bin/env bash 
set -e

echo "pyPamtra2 installing script"
echo "USAGE: ./install.sh [-p|--python PYTHON_COMMAND] [-f|--flag FLAG1 -f|--flag FLAG2 -f|--flag FLAG3 ...]"
echo "TYPICAL command: ./install.sh -p python3 -f --user"
echo "This will install the pyPamtra2 subpackages under /$HOME/.local/lib/python3.X/ which is already in your PYTHONPATH on most unix systems"
echo "This script is not robust enough to accept any parameter. Anything different from -p -f cause the script to freeze"

PYTHON='python'

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -p|--python)
    PYTHON="$2"
    shift # past argument
    shift # past value
    ;;&
    -f|--flag)
    FLAGS+=" $2"
    shift # past argument
    shift # past value
    ;;&
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

echo "python command " $PYTHON
echo "passed flags " $FLAGS

cd pyPamtraRadarSimulator
$PYTHON setup.py install $FLAGS
cd ../

cd pyPamtraRadarMoments
$PYTHON setup.py install $FLAGS
cd ../

cd pyPamtra2
$PYTHON setup.py install $FLAGS
cd ../

cd refractive
$PYTHON setup.py install $FLAGS
cd ../

cd scattering
$PYTHON setup.py install $FLAGS
cd ../
