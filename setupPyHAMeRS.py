#!/usr/bin/env python

from numpy.distutils.core import setup, Extension

import sys
import os
import subprocess

# Python version
if sys.version_info[:2] < (3, 6):
    print('PyHAMeRS requires Python 3.6 or newer')
    sys.exit(-1)

version = '1.0'

# Modules.
modules = ['pyhamers.derivatives',
           'pyhamers.derivatives.explicit',
           'pyhamers.upsampling',
           'pyhamers.readers']

setup(name = 'PyHAMeRS',
      version = version,
      description = 'Postprocessing utilities for HAMeRS',
      author = 'Contributors of HAMeRS',
      author_email = 'wongml@stanford.edu',
      packages = ['pyhamers'] + modules,
      )