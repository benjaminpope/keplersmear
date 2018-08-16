from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import Configuration
import distutils.sysconfig as ds

long_description = "This package generates light curves from collateral 'smear' data recorded by Kepler and K2. \
This is useful for reconstructing light curves of stars that were otherwise not observed by the Kepler mission because they were too bright."

setup(name='keplersmear',
      version='0.5',
      description='Kepler and K2 Smear Photometry.',
      long_description=long_description,
      author='Benjamin Pope',
      author_email='benjamin.pope@physics.ox.ac.uk',
      url='',
      package_dir={'keplersmear':'src'},
      scripts=['bin/kepler_smear_list','bin/kepler_smear_lc','bin/k2_smear_lc'],
      packages=['keplersmear'],
      install_requires=["numpy", "astropy", "scipy"],
      license='GPLv3',
      classifiers=[
          "Topic :: Scientific/Engineering",
          "Intended Audience :: Science/Research",
          "Intended Audience :: Developers",
          "Development Status :: 4 - Beta",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
          "Operating System :: OS Independent",
          "Programming Language :: Python"
      ]
     )
