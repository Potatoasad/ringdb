#!/usr/bin/env python

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

description = '''A gravitational event database that queries and locally saves event strain samples, 
detector PSDs and posterior samples - all in one place'''

setup(name='ringdb',
      version='0.1',
      description=description.replace('\n',''),
      author='Asad Hussain',
      author_email='asadh@utexas.edu',
      url='https://github.com/potatoasad/ringdb',
      download_url='https://github.com/Potatoasad/ringdb/archive/refs/tags/v0.1.0.zip',
      license='MIT',
      packages=['ringdb'],
      package_data={'ringdb': ['metadb/*']},
      install_requires=[
            'Cython>=0.22',
            'arviz',
            'h5py',
            'lalsuite',
            'matplotlib',
            'numpy',
            'pandas',
            'pystan>=2,<3',
            'qnm',
            'scipy',
            'seaborn',
            'gwosc',
            'ringdown']
     )
