#!/usr/bin/python

from distutils.core import setup
from castepconv import __version__

setup(name='castepconv',
    version=__version__,
    packages = ['cconv'],
    scripts=['castepconv.py'],
    )
