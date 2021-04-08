#!/usr/bin/env python

from distutils.core import setup

LONG_DESCRIPTION = \
'''Find processed pseudo genes in DNA sequencing data using input
structural variant calls'''


setup(
    name='pseudofinder',
    version='0.1.0.0',
    author='Bernie Pope',
    author_email='bjpope@unimelb.edu.au',
    packages=['pseudofinder'],
    package_dir={'pseudofinder': 'pseudofinder'},
    entry_points={
        'console_scripts': ['pseudofinder = pseudofinder.pseudofinder:main']
    },
    url='https://github.com/GITHUB_USERNAME/pseudofinder',
    license='LICENSE',
    description=('Find processed pseudo genes in DNA sequencing data using input structural variant calls'),
    long_description=(LONG_DESCRIPTION),
    install_requires=["intervaltree", "cyvcf2"],
)
