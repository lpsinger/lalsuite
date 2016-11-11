from setuptools import setup

setup(
  name = "glue",
  version = "1.54.0",
  author = "Leo Singer",
  author_email = "leo.singer@ligo.org",
  description = "Transitional stub for ligo-glue",
  url = "http://www.lsc-group.phys.uwm.edu/daswg/",
  license = 'GPL-v3',
  packages = ['glue'],
  install_requires = ['ligo-glue > 1.53.0'],
  classifiers = [
    'Development Status :: 5 - Production/Stable',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Operating System :: POSIX',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3.4',
    'Topic :: Scientific/Engineering :: Astronomy',
    'Topic :: Scientific/Engineering :: Physics'
  ]
)
