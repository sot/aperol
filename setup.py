# Licensed under a 3-clause BSD style license - see LICENSE.rst
from setuptools import setup

entry_points = {'console_scripts': ['aperoll=aperoll.scripts.aperoll:main']}

setup(name='aperoll',
      author='Javier Gonzalez',
      description='Prosecco with Aperoll make a nice spritz',
      author_email='javier.gonzalez@cfa.harvard.edu',
      packages=['aperoll', 'aperoll.scripts', 'aperoll.widgets'],
      license=("New BSD/3-clause BSD License\nCopyright (c) 2021"
               " Smithsonian Astrophysical Observatory\nAll rights reserved."),
      url='http://www.github.com/sot/aperoll',
      entry_points=entry_points,
      zip_safe=False,
      use_scm_version=True,
      setup_requires=['setuptools_scm', 'setuptools_scm_git_archive'],
      )
