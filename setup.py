#!/usr/bin/env python

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

description = '''A gravitational event database that queries and locally saves event strain samples, 
detector PSDs and posterior samples - all in one place'''

setup(name='ringdb',
      version='0.1.2',
      description=description.replace('\n',''),
      author='Asad Hussain',
      author_email='asadh@utexas.edu',
      url='https://github.com/potatoasad/ringdb',
      download_url='https://github.com/Potatoasad/ringdb/archive/refs/tags/v0.1.2.zip',
      license='MIT',
      packages=['ringdb'],
      package_data={'ringdb': ['metadb/*']},
      install_requires=[
            'h5py',
            'lalsuite',
            'matplotlib',
            'numpy',
            'pandas',
            'ringdown']
     )
