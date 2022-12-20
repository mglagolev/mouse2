#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

    Calculate the molecular ordering parameters for lamellae
    containing tilted copolymer blocks, as described in the
    paper by M. A. Osipov, M. V. Gorkunov, A. V. Berezkin,
    A. A. Antonov and Y. V. Kudryavtsev "Molecular theory
    of the tilting transition and computer simulations of
    the tilted lamellar phase of rod–coil diblock copolymers"
    https://doi.org/10.1063/5.0005854
    

@author: Anna Glagoleva, Mikhail Glagolev

"""

import MDAnalysis as mda
import argparse
import json
from lib.lamellar_orientation import lamellar_orientational_ordering_parameters


if __name__ == "__main__":
    

    parser = argparse.ArgumentParser(
        description = 'Calculate the lamellar orientational ordering'
        + ' parameters')

    parser.add_argument(
        'input', metavar = 'INPUT', action = "store", nargs = '+',
        help = "input files")

    parser.add_argument(
        '--block-types', metavar = 'TYPES', type = str, nargs = 2,
        default = ["1", "2"],
        help = "bead types for the blocks A dnd B (provide 0 or 2 arguments)")
    
    parser.add_argument(
        '--A', action = "store_true", help = "Calculate the values for "
                                           + "block A")
    
    parser.add_argument(
        '--B', action = "store_true", help = "Calculate the values for "
                                           + "block B")
    
    parser.add_argument(
        '--verbose', action = "store_true", help = "Store the values for "
                                                 + "individual molecules")
    

    args = parser.parse_args()

    u = mda.Universe(*args.input)
    
    result = lamellar_orientational_ordering_parameters(u, 
                                             args.block_types[0],
                                             args.block_types[1],
                                             store_A_values = args.A,
                                             store_B_values = args.B,
                                             store_block_values = args.verbose
                                             )
    print(json.dumps(result, indent = 2))