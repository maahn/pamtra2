import os
import sys
from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def configuration(parent_package='', top_path=None):

    config = Configuration(
        'pyPamtraRadarSimulator', parent_package, top_path,
        version='0.1',
        author="Maximilin Maahn",
        author_email="maximilian.maahn@colorado.edu",
        description="Simulate a radar Doppler spectrum",
        license="GPL v3",
        python_requires='>=3.5',
        # url='https://github.com/maahn/pyPamtraRadarSimulator',
        # download_url='https://github.com/maahn/pyPamtraRadarSimulator/releases/download/0.1/pyPamtraRadarSimulator-0.1.zip',
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
    config.add_extension('pyPamtraRadarSimulatorLib',
                         sources=[
                             'pyPamtraRadarSimulator/pyPamtraRadarSimulatorLib/pyPamtraRadarSimulatorLib.pyf',
                             'pyPamtraRadarSimulator/pyPamtraRadarSimulatorLib/kinds.f90',
                             'pyPamtraRadarSimulator/pyPamtraRadarSimulatorLib/report_module.f90',
                             'pyPamtraRadarSimulator/pyPamtraRadarSimulatorLib/constants.f90',
                             'pyPamtraRadarSimulator/pyPamtraRadarSimulatorLib/dsort.f90',
                             'pyPamtraRadarSimulator/pyPamtraRadarSimulatorLib/convolution.f90',
                             'pyPamtraRadarSimulator/pyPamtraRadarSimulatorLib/viscosity_air.f90',
                             'pyPamtraRadarSimulator/pyPamtraRadarSimulatorLib/random_module.f90',
                             'pyPamtraRadarSimulator/pyPamtraRadarSimulatorLib/rescale_spectra.f90',
                             'pyPamtraRadarSimulator/pyPamtraRadarSimulatorLib/dia2vel.f90',
                             'pyPamtraRadarSimulator/pyPamtraRadarSimulatorLib/radar_simulator.f90',
                             'pyPamtraRadarSimulator/pyPamtraRadarSimulatorLib/radar_spectral_broadening.f90',
                             'pyPamtraRadarSimulator/pyPamtraRadarSimulatorLib/radar_spectrum.f90',
                             'pyPamtraRadarSimulator/pyPamtraRadarSimulatorLib/rho_air.f90',
                         ],
                         library_dirs=['~/.local/lib/', '/usr/local/lib/'],
                         libraries=['fftw3', 'lapack'],
                         **kw)

    return config


if __name__ == "__main__":

    setup(configuration=configuration,
          packages=['pyPamtraRadarSimulator',
                    'pyPamtraRadarSimulator.pyPamtraRadarSimulatorLib'],
          # package_data = {
          #     'pyPamtraRadarSimulator': ['file'],
          # },
          platforms=['any'],
          requires=['numpy', 'scipy'])
