import os
import sys
from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup




def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def configuration(parent_package='',top_path=None):
    
    config = Configuration('pyPamtraRadarMoments', parent_package, top_path,
        version = '0.1',
        author  = "Maximilin Maahn",
        author_email = "maximilian.maahn@colorado.edu",
        description = "Estimate Moments from a radar Doppler spectrum",
        license = "GPL v3",
        python_requires='>=3.5',
        url = 'https://github.com/maahn/pyPamtraRadarMoments',
        download_url = 'https://github.com/maahn/pyPamtraRadarMoments/releases/download/0.1/pyPamtraRadarMoments-0.1.zip',
        long_description = read('README.rst'),
        classifiers = [
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
    config.add_extension('pyPamtraRadarMomentsLib',
        sources=[
                 'pyPamtraRadarMoments/pyPamtraRadarMomentsLib/nan.f90',
                 'pyPamtraRadarMoments/pyPamtraRadarMomentsLib/pyPamtraRadarMomentsLib.pyf',
                 'pyPamtraRadarMoments/pyPamtraRadarMomentsLib/kinds.f90',
                 'pyPamtraRadarMoments/pyPamtraRadarMomentsLib/report_module.f90',
                 'pyPamtraRadarMoments/pyPamtraRadarMomentsLib/constants.f90',
                 'pyPamtraRadarMoments/pyPamtraRadarMomentsLib/convolution.f90',
                 'pyPamtraRadarMoments/pyPamtraRadarMomentsLib/dsort.f90',
                 'pyPamtraRadarMoments/pyPamtraRadarMomentsLib/hildebrand_sekhon.f90',
                 'pyPamtraRadarMoments/pyPamtraRadarMomentsLib/smooth_savitzky_golay.f90',
                 'pyPamtraRadarMoments/pyPamtraRadarMomentsLib/calc_moments.f90',
                 ],
        library_dirs = ['/usr/local/lib/'],
        libraries = ['fftw3','lapack'],
        **kw)

    return config


if __name__ == "__main__":

    
    setup(configuration=configuration,
        packages = ['pyPamtraRadarMoments','pyPamtraRadarMoments.pyPamtraRadarMomentsLib'],        
        # package_data = {
        #     'pyPamtraRadarMoments': ['file'],
        # },
        platforms = ['any'],
        requires = ['numpy', 'scipy'])

