# -*- coding: utf-8 -*-
"""
    Copyright (C) 2017 - 2018 Davide Ori dori@uni-koeln.de
    Institute for Geophysics and Meteorology - University of Cologne

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import os
import sys
import numpy

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

from distutils.core import setup

# Build the main package, with script etc...
# ==========================================
setup(name = 'singleScattering',
      description = 'Set of utilities to compute the scattering properties of particles',
      author      = 'Davide Ori',
      author_email = 'dori@uni-koeln.de',
      packages    = ['singleScattering'],
      license = 'GPL v3', 
      python_requires='>=3.5', 
      #scripts = ["scripts/myscript.py"],
      version = '0.1',
      url = 'https://github.com/maahn/pyPamtra2',
      download_url = 'https://github.com/maahn/pyPamtra2/releases/download/0.1/pyPamtra2-0.1.zip',
      long_description = read('README.rst'),
      classifiers = [
            "Development Status :: 3 - Alpha",
            "License :: OSI Approved :: GPL v3 License",
            "Operating System :: OS Independent",
            "Programming Language :: Python",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering :: Atmospheric Science",
            ]
      )


# Build the extensions
# --------------------

# Build the cython extension
from distutils.extension import Extension
from Cython.Distutils import build_ext

cMie = Extension("singleScattering.cMie",
                sources=["Mie/cython/cMie.pyx",
                         "Mie/src/cMie.c"],
                include_dirs=[numpy.get_include()],
                extra_compile_args=["-O3", "-ffast-math", "-Wall", "-lm", "-fPIC", "-std=c99"],
                language='c')

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [cMie]
    )

# Build the fortran modules with f2py
"""
# Build the f2py fortran extension
# --------------------------------
from numpy.distutils.core import Extension
from numpy.distutils.core import setup

fortranlib = Extension(name = 'mypackage.flib', # you will call in python mypackage.flib.FORTRAN_IMPLEMENTED_FUNCTION()
                 extra_compile_args = ['-O3'],
                 sources = ['src_fortran/mymodule.f90'], # you may add several modules files under the same extension
                 )

setup(
    ext_modules = [fortranlib]
)
"""
