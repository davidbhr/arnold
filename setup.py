#!/usr/bin/env python

from setuptools import setup

setup(
    name='arnold',
    packages=['arnold'],
    version='1.0',
    description='Python package to analyze musccle strength using 3D Traction Force Microscopy',
    url='',
    download_url = '',
    author='David Böhringer, Christoph Mark',
    author_email='davidboehringe@gmail.com',
    license='The MIT License (MIT)',
    install_requires=['numpy>=1.16.2',
                      'pandas>=0.23.4',
                      'matplotlib>=2.2.2',
					  'roipoly>=0.5.2'],
    keywords = ['muscle', 'contractility', 'material simulation', 'Traction Force Microscopy'],
    classifiers = [],
    )
