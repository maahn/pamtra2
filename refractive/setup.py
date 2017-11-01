import os
import sys
from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def configuration(parent_package='',top_path=None):
    
    config = Configuration('pyPamtra2', parent_package, top_path,
        version = '0.1',
        author  = "Pamtra Team",
        author_email = "dori@uni-koeln.de",
        description = "complex refractive index of ice and water",
        license = "GPL v3",
        python_requires='>=3.5',
        url = 'https://github.com/maahn/pyPamtra2',
        download_url = 'https://github.com/maahn/pyPamtra2/releases/download/0.1/pyPamtra2-0.1.zip',
        long_description = read('README.rst'),
        classifiers = [
            "Development Status :: 3 - Alpha",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            #"Programming Language :: Fortran",
            "Programming Language :: Python",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering :: Atmospheric Science",
            ]
    )

    return config


if __name__ == "__main__":
    
    setup(configuration=configuration,
        packages = ['pyPamtra2'],        
        # package_data = {
        #     'pyPamtra2': ['file'],
        # },
        platforms = ['any'],
        requires = ['numpy', 'scipy','xarray','json'])

