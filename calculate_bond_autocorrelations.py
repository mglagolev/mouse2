#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 19:59:20 2022

@author: Mikhail Glagolev
"""
import MDAnalysis as mda
from MDAnalysis import transformations
import numpy as np

def calculate_bond_autocorrelations(u: mda.Universe, k_max,
                                    selection = None,
                                    different_molecules: bool = False):
    """
    
    Calculate the autocorrelation function of the polymer bonds.
    The formula is presented in https://doi.org/10.1134/S0965545X10070102
    Application in all-atom simulations: https://doi.org/10.3390/polym11122056

    Parameters
    ----------
    u : mda.Universe
        DESCRIPTION.
    k_max : integer
        The maximum value of the distance between the bonds along the backbone
    different_molecules : bool
        Take into account the bonds where the particles have different
        residue ids. The default value is False.

    Returns
    -------
    For every timestep a list of float values [C(0), C(1), ...]
    The lists are wrapped into a list for all the timesteps.

    """
    # Unwrap all the coordinates, so that all the bond lengths are correct.
    unwrap = transformations.unwrap(u.atoms)
    u.trajectory.add_transformations(unwrap)
    
    # Select atoms by type or read selection criteria in MDAnalysis synthax
    if selection is not None:
        atoms = u.select_atoms(selection)
    else:
        atoms = u.atoms
    
    results = []
    
    for ts in u.trajectory:
        # List of c[k] = c(k), k = 0, 1, 2, ...
        ck = []
        # Total number of bonds
        nbonds = len(atoms.bonds)
        # Molecule ids
        bond_resids = atoms.bonds.atom1.resids
        #Determine the vectors for all the bonds with NumPy
        b = atoms.bonds.atom2.positions - atoms.bonds.atom1.positions   
        # Consider each value of the shift between the bonds k
        for k in range(0, k_max + 1):
            # Create two array shifted by k, by padding the arrays with
            # the values of "1." (so that vector length is not zero)
            b1 = np.concatenate((b, np.full((k, 3), 1.)), axis = 0)
            b2 = np.concatenate((np.full((k, 3), 1.), b), axis = 0)
            # Consider the calculation result valid if neither value is padded
            valid1 = np.concatenate((np.full((nbonds,), True),
                                    np.full((k,), False)))
            valid2 = np.concatenate((np.full((k,), False),
                                    np.full((nbonds,), True)))
            valid = np.logical_and(valid1, valid2)
            # If cross-molecule correlations should not be accounted for,
            # then also check for the equality of the resids
            if not different_molecules:
                # Pad the residue id arrays
                resid1 = np.concatenate((bond_resids, np.zeros((k,))))
                resid2 = np.concatenate((np.zeros((k,)), bond_resids))
                # Take into account only molecules with same residue id
                valid = np.logical_and(valid, np.equal(resid1, resid2))
            # Mask is True for the values that are not valid    
            mask = np.logical_not(valid)
            # Calculate the correlation values for all the bonds
            c = ( np.sum( np.multiply(b1, b2), axis = 1)
                 / np.linalg.norm(b1, axis = 1)
                 / np.linalg.norm(b2, axis = 1))
            c_masked = np.ma.masked_array(c, mask = mask)
            c_average = np.ma.average(c_masked)
            ck.append(c_average)
        results.append(ck)
    return results

if __name__ == "__main__":
    
    import argparse

    parser = argparse.ArgumentParser(
        description = 'Calculate backbone bonds autocorrelations')

    parser.add_argument(
        'input', metavar = 'INPUT', action = "store", nargs = '+',
        help = "input files")

    parser.add_argument(
        '--k_max', metavar = 'k_max', type = int, nargs = '?',
        default = 0.,
        help = "maximum distance between the bonds along the backbone")
    
    parser.add_argument(
        '--selection', metavar = 'QUERY', type = str, nargs = '?',
        help 
        = "Consider only selected atoms, use MDAnalysis selection language")
    
    parser.add_argument(
        "--different-molecules", action = "store_true",
        help = "Calculate corrlations based on particle index number,\
            even if the bonds belong to different molecules")
    

    args = parser.parse_args()

    u = mda.Universe(*args.input)
    
    ck_series = calculate_bond_autocorrelations(u, k_max = args.k_max,
                                                selection = args.selection,
                                                different_molecules =
                                                args.different_molecules)
    
    # Print the values:
    for ck_values in ck_series:
        for k in range(len(ck_values)):
            print(str(k) + " " + str(ck_values[k]))
        print("\n")