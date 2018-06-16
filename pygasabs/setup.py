import os
import sys
from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def configuration(parent_package='', top_path=None):

    config = Configuration(
        'pygasabs', parent_package, top_path,
        version='0.1',
        author="Maximilin Maahn",
        author_email="maximilian.maahn@colorado.edu",
        description="Estimates gaseous attenuation in the microwave region.",
        license="GPL v3",
        python_requires='>=3.5',
        url='https://github.com/maahn/pamtra2/pygasabs',
        # download_url='https://github.com/maahn/pygasabs/releases/download/0.1/pygasabs-0.1.zip',
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

    kw = {}
    if sys.platform == 'darwin':
        kw['extra_link_args'] = ['-undefined dynamic_lookup', '-bundle']
    config.add_extension('pygasabs_lib',
                         sources=[
                             'pygasabs/pygasabs_lib/pygasabs_lib.pyf',
                             'pygasabs/pygasabs_lib/kinds.f90',
                             'pygasabs/pygasabs_lib/report_module.f90',
                             'pygasabs/pygasabs_lib/constants.f90',
                             'pygasabs/pygasabs_lib/gasabs_module.f90',
                             'pygasabs/pygasabs_lib/rosen98_gasabs.f90',
                             'pygasabs/pygasabs_lib/mpm93.f90',
                         ],
                         # library_dirs=['/usr/local/lib/'],
                         # libraries=['fftw3', 'lapack'],
                         **kw)

    return config

if __name__ == "__main__":

    setup(configuration=configuration,
          packages=['pygasabs', 'pygasabs.pygasabs_lib'],
          # package_data = {
          #     'pygasabs': ['file'],
          # },
          platforms=['any'],
          requires=['numpy', 'scipy'])
