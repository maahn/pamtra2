import os
import shlex

import sys
from glob import glob

import numpy
import setuptools
from Cython.Build import cythonize
from numpy.distutils.core import setup
from numpy.distutils.extension import Extension
from numpy.distutils.misc_util import Configuration


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


kw = {}
if sys.platform == 'darwin':
    kw['extra_link_args'] = ['-undefined dynamic_lookup', '-bundle']


def configuration(parent_package='', top_path=None):

    config = Configuration(
        'pamtra2',
        parent_package,
        top_path,
        version='0.1',
        author="Pamtra Team",
        author_email="maximilian.maahn@colorado.edu",
        description="atmospheric microwace passive"
                    " and active instrument simulator",
        license="GPL v3",
        python_requires='>=3.5',
        url='https://github.com/maahn/pamtra2',
        download_url='https://github.com/maahn/pamtra2/'
                     'releases/download/0.1/pamtra2-0.1.zip',
        long_description=read('README.rst'),
        classifiers=[
            "Development Status :: 3 - Alpha",
            "License :: OSI Approved :: GPL v3 License",
            "Operating System :: OS Independent",
            "Programming Language :: Fortran",
            "Programming Language :: Python",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering :: Atmospheric Science",
        ]
    )

    return config


refractiveIndex_path = 'libs/refractiveIndex/refractiveIndex'
meteo_si_path = 'libs/meteo_si/meteo_si'

singleScattering_path = 'libs/singleScattering'
cMie = Extension(
    name = "pamtra2.libs.singleScattering.cMie",
    sources = ["%s/Mie/cython/cMie.pyx" % singleScattering_path,
             "%s/Mie/src/cMie.c" % singleScattering_path],
    include_dirs = [numpy.get_include()],
    extra_compile_args = ["-O3", "-ffast-math",
                          "-Wall", "-lm", "-fPIC", "-std=c99"],
    language='c'
)


if __name__ == "__main__":

    setup(
        configuration=configuration,
        packages=['pamtra2', 'pamtra2.hydrometeors',
                  'pamtra2.instruments', 'pamtra2.libs',
                  'pamtra2.libs.refractiveIndex',
                  'pamtra2.libs.singleScattering',
                  'pamtra2.libs.meteo_si','pamtra2.importer'
                  ],
        package_dir={
            'pamtra2': 'src/pamtra2' ,
            # 'pamtra2.libs.pamgasabs': pamgasabs_path,
            # 'pamtra2.libs.pyPamtraRadarSimulator': pyrasim_path,
            # 'pamtra2.libs.pyPamtraRadarMoments': pyramom_path,
            'pamtra2.libs.refractiveIndex': refractiveIndex_path,
            'pamtra2.libs.meteo_si': meteo_si_path,
            'pamtra2.libs.singleScattering': '%s/singleScattering' %
            singleScattering_path,
        },
        package_data={
            # Don't aks me why, but I need an extra random character in front
            # of the *dat files...
            'pamtra2.libs.refractiveIndex': ['*.dat'],
            'pamtra2.libs.singleScattering': ['*.dat'],
        },
        platforms=['any'],
        install_requires=['numpy', 'scipy', 'xarray', 'dask[complete]', 'numba',
                  'cython'],
        build_requires=['numpy', 'cython'],
        setup_requires=["pytest-runner"],
        tests_require=["pytest"],
        ext_modules=cythonize(
            [cMie]),

    )
