#!/usr/bin/env python

from setuptools import setup

setup(
    name='radpipe',
    version='0.0.1',
    author='Jessica Chung, Bernie Pope',
    author_email='jchung@unimelb.edu.au',
    packages=['radpipe', 'radpipe.pipeline_base'],
    entry_points={
        'console_scripts': ['radpipe = radpipe.main:main']
    },
    url='https://github.com/jessicachung/radpipe',
    license='LICENSE',
    description='A pipeline for RAD-seq analysis',
    long_description=open('README.md').read(),
    install_requires=[
        "ruffus >= 2.6.3",
        "drmaa == 0.7.6",
        "PyYAML == 3.11"
    ],
)
