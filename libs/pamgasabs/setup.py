import os
import sys

from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def configuration(parent_package='', top_path=None):

    config = Configuration(
        'pamgasabs', parent_package, top_path,
        version='0.1',
        author="Maximilin Maahn",
        author_email="maximilian.maahn@colorado.edu",
        description="Estimates gaseous attenuation in the microwave region.",
        license="GPL v3",
        python_requires='>=3.5',
        url='https://github.com/maahn/pamtra2/pamgasabs',
        # download_url='https://github.com/maahn/pamgasabs/releases/download/0.1/pamgasabs-0.1.zip',
        long_description=read('README.rst'),
        classifiers=[
           "Development Status :: 3 - Alpha",
           "License :: OSI Approved :: MIT License",
           "Operating System :: OS Independent",
           "Programming Language :: Fortran",
           "Programming Language :: Python",
           "Intended Audience :: Science/Research",
           "Topic :: Scientific/Engineering :: Atmospheric Science",
        ]
        )

    kw = {
    'extra_compile_args' : '-Wall -fcheck=all -O3'
    }

    if sys.platform == 'darwin':
        kw['extra_link_args'] = ['-undefined dynamic_lookup', '-bundle']
    config.add_extension('pamgasabs_lib',
                         sources=[
                             'pamgasabs/pamgasabs_lib/pamgasabs_lib.pyf',
                             'pamgasabs/pamgasabs_lib/kinds.f90',
                             'pamgasabs/pamgasabs_lib/report_module.f90',
                             'pamgasabs/pamgasabs_lib/constants.f90',
                             'pamgasabs/pamgasabs_lib/gasabs_module.f90',
                             'pamgasabs/pamgasabs_lib/rosen98_gasabs.f90',
                             'pamgasabs/pamgasabs_lib/mpm93.f90',
                         ],
                         # library_dirs=['/usr/local/lib/'],
                         # libraries=['fftw3', 'lapack'],
                         **kw)

    return config

if __name__ == "__main__":

    setup(configuration=configuration,
          packages=['pamgasabs', 'pamgasabs.pamgasabs_lib'],
          # package_data = {
          #     'pamgasabs': ['file'],
          # },
          platforms=['any'],
          requires=['numpy', 'scipy'])
