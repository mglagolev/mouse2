#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Mikhail Glagolev
"""
from setuptools import setup

from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='mouse2',
    version='0.1.8',    
    description="""A toolkit for processing molecular dynamics simulation data
    with a focus on chiral ordering""",
    url='https://github.com/mglagolev/mouse2',
    author='Mikhail Glagolev, Anna Glagoleva',
    author_email='mikhail.glagolev@gmail.com',
    license='GNU GPL v3',
    packages=['mouse2', 'mouse2.lib'],
    install_requires=['numpy',
                      'MDAnalysis',
                      'networkx',
                      ],
    entry_points = {'console_scripts': ['aggregates = mouse2.aggregates:main',
    'backbone_bond_autocorrelations = mouse2.backbone_bond_autocorrelations:main',
    'backbone_twist = mouse2.backbone_twist:main',
    'bond_orientational_ordering = mouse2.bond_orientational_ordering:main',
    'lamellar_orientational_ordering = mouse2.lamellar_orientational_ordering:main',
    'data2pdb = mouse2.data2pdb:main',
    ]},
    long_description = long_description,
    long_description_content_type='text/markdown',

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',  
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 3',
    ],
)
