#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Mikhail Glagolev
"""
from setuptools import setup

setup(
    name='mouse2',
    version='0.1.1',    
    description="""A toolkit for processing molecular dynamics simulation data
    with a focus on chiral ordering""",
    url='https://github.com/mglagolev/mouse2',
    author='Mikhail Glagolev, Anna Glagoleva',
    author_email='mikhail.glagolev@gmail.com',
    license='GNU GPL v3',
    packages=['mouse2'],
    install_requires=['numpy',
                      'MDAnalysis',
                      'networkx',
                      ],

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',  
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 3',
    ],
)
