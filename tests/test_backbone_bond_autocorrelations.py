#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 19:50:56 2023

@author: Mikhail Glagolev
"""

import unittest
import MDAnalysis as mda
import os
from mouse2.backbone_bond_autocorrelations import calculate_bond_autocorrelations
import numpy as np

test_dir = absolute_path = os.path.dirname(__file__)

test_file = os.path.join(test_dir, "helical_lamellae.data.gz")

tolerance = 1e-6

target_flexible_ck = np.asarray([
     1.0,
     0.204648203125,
     0.08333987862723215,
     0.04822373422475962,
     0.04239111735026042,
     0.028410007546164774,
     0.0326979833984375,
     0.0241927001953125,
     0.019064814758300783,
     0.017419686453683036,
     0.017799589029947917,
])

target_helical_ck = np.asarray([
     1.0,
     -0.016461666666666666,
     -0.2911587890625,
     0.8378388221153846,
     0.3077317708333333,
     -0.41193316761363635,
     0.607004609375,
     0.5848464409722223,
     -0.376908251953125,
     0.31830694754464284,
     0.747251171875
])

class TestAutocorrelations(unittest.TestCase):

    def test_flexible_ck(self):
        u = mda.Universe(test_file)
        result = calculate_bond_autocorrelations(u, 10, selection = "type 1")
        ck = np.asarray(list(result["data"].values())[0])
        discrepancy = np.max(np.abs(ck - target_flexible_ck))
        assert discrepancy <= tolerance
    
    def test_helical_ck(self):
        u = mda.Universe(test_file)
        result = calculate_bond_autocorrelations(u, 10, selection = "type 2")
        ck = np.asarray(list(result["data"].values())[0])
        discrepancy = np.max(np.abs(ck - target_helical_ck))
        assert discrepancy <= tolerance
        
if __name__ == "__main__":
    unittest.main()