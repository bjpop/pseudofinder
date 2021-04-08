#!/usr/bin/env python

from distutils.core import setup

LONG_DESCRIPTION = \
'''Find processed pseudo genes in DNA sequencing data using input
structural variant calls'''


setup(
    name='psuedofinder',
    version='0.1.0.0',
    author='Bernie Pope',
    author_email='bjpope@unimelb.edu.au',
    packages=['psuedofinder'],
    package_dir={'psuedofinder': 'psuedofinder'},
    entry_points={
        'console_scripts': ['psuedofinder = psuedofinder.psuedofinder:main']
    },
    url='https://github.com/GITHUB_USERNAME/psuedofinder',
    license='LICENSE',
    description=('Find processed pseudo genes in DNA sequencing data using input structural variant calls'),
    long_description=(LONG_DESCRIPTION),
    install_requires=["intervaltree"],
)
