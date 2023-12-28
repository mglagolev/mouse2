#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a sample script to demonstrate the usage of the MOUSE2
local_alignment module.

To run the script, the MOUSE2 package shall be installed either using pip:
    pip install mouse2
or by downloading from the GitHub repository:
    https://github.com/mglagolev/mouse2

In the latter case, the following dependencies shall be installed by the user
manually:
    - MDAnalysis
    - numpy
    - networkx
    - matplotlib
    - scipy
and the import command needs to be modified accordingly.

The identifiers for the atoms at the ends of the helical fragments were in
this case identified by us manually and are stored in the atom_nums variable.

The script is run without any arguments, the lower and upper cutoff radii
for the analysis can be modified through the r_mins and r_maxes variables.
"""

from urllib.request import urlretrieve
import MDAnalysis as mda
import numpy as np
import tempfile
import os
from mouse2.local_alignment import local_alignment

atom_nums = [
    [5,12],
    [14,29],
    [34,54],
    [64,74],
    [83,91],
    [94,103],
    [119,127],
    [129,144],
    [150,167],
    [174,204],
    [206,221],
    [225,245],
    [248,266],
    [323,336],
    [341,360],
    [364,395],
    [398,413],
    [418,437],
    [443,464],
    [516,535],
    [539,559],
    [562,582],
]


atom_ids = [[i[0]-1, i[1]-1] for i in atom_nums]

r_maxes = np.arange(8, 32, dtype = int)
r_mins = np.repeat(1e-6, len(r_maxes))
#r_mins = np.arange(1, 13, dtype = int)

_, tmpfilename = tempfile.mkstemp(suffix = '.pdb')

#Set the maximum correlation length

#Retreive the PDB file
url = "https://files.rcsb.org/download/6OCK.pdb"
filename = "6OCK.pdb"
urlretrieve(url, filename)

u = mda.Universe(filename)

b = u.select_atoms('backbone and name N')
ow = b.split('segment')
bonds = []
bond_types = []
for i in ow:
    for j in range(len(i)-1):
        bonds.append([i[j].ix, i[j+1].ix])
        bond_types.append('1')
u.add_bonds(bonds, types = bond_types)
b = u.select_atoms('backbone and name N')
b.write(tmpfilename)

print("# r\ts")
for i in range(len(r_mins)):
    r_min, r_max = r_mins[i], r_maxes[i]
    #r = (r_min + r_max) / 2.
    r = r_max
    u = mda.Universe(tmpfilename)
    result = local_alignment(u, r_min = r_min, r_max = r_max,
                             mode = 'average', id_pairs = atom_ids)
    print(f"{r:.1f}\t{list(result['data'].values())[0]['average_s']:.3f}")

os.remove(tmpfilename)