#!/usr/bin/python

from castepconv import __version__
from setuptools import setup, find_packages

if __name__ == "__main__":

    setup(name='CastepConv',
          version=__version__,
          description='A tool to automate convergence with CASTEP',
          url='https://github.com/CCP-NC/castepconv',
          author='Simone Sturniolo',
          author_email='simone.sturniolo@stfc.ac.uk',
          license='GPL',
          packages=['cconv'],
          scripts=['castepconv.py'],
          install_requires=[
              'numpy',
              'ase>=3.17'
          ],
          python_requires='>=2.6'
          )
