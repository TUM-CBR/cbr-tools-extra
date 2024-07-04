#!/usr/bin/env python

from distutils.core import setup

setup(
    name='cbr-tools-extra',
    version='0.1.0',
    description='Support module for cbr-tools.',
    author='netogallo',
    author_email='connect@netowork.me',
    url='https://github.com/TUM-CBR/cbr-tools-extra',
    packages=['cbrextra'],
    scripts=['./cbrtools'],
    install_requires = [
       'biopython',
       'docopt',
       'igraph',
       'openpyxl',
       'open3d-cpu',
       'pandas',
       'primer3-py',
       'pydantic',
       'requests',
       'SQLAlchemy'
    ],
)
