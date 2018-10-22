#!/usr/bin/python

from distutils.core import setup
from castepconv import __vers_number__

setup(name='castepconv',
    version=__vers_number__,
    packages = ['cconv'],
    scripts=['castepconv.py'],
    )
