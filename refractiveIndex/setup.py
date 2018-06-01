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
from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def configuration(parent_package='',top_path=None):
    
    config = Configuration('refractiveIndex', parent_package, top_path,
        version = '0.1',
        author  = "Davide Ori",
        author_email = "dori@uni-koeln.de",
        description = "complex refractive index of ice and water",
        license = "GPL v3",
        python_requires='>=3.5',
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

    return config


if __name__ == "__main__":
    
    setup(configuration=configuration,
        packages = ['refractiveIndex'],        
        package_data = {
            'refractiveIndex': ['*.dat'],
        },
        platforms = ['any'],
        requires = ['numpy','scipy'])

