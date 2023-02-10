#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 17:20:41 2022

@author: Mikhail Glagolev
"""

import MDAnalysis as mda
from MDAnalysis import transformations
import numpy as np
from mouse2.lib.vector_orientational_ordering \
    import calculate_cos_sq_for_reference
import json

def averaged_frequencies_bin_centers(result, frequencies_key, bin_edges_key):
    """
    Return averaged histogram data across the time steps.
    Bin edges are assumed to be the same for all the timesteps.

    """
    bin_edges = np.asarray(list(result["data"].values())[0][bin_edges_key])
    bincenters = (bin_edges[1:] + bin_edges[:-1]) / 2.
    n_bins = len(bincenters)
    frequencies = np.ndarray((n_bins))
    for ts in result["data"]:
        frequencies += np.asarray(
                    result["data"][ts][frequencies_key])
    frequencies /= len(result["data"])
    return frequencies, bincenters
    

def calculate_orientation_order_parameter(
        u: mda.Universe, r_min = 0., r_max = -1., n_bins = 150,
        mode = 'average', same_molecule = True):
    """
    
    This function calculates the angles between the bonds, if their
    midpoints are located within the range of [rmin, rmax].
    The local ordering parameter is then calculated as
    S = 3/2 ( <(cos(gamma))^2>) - 1/2)
    where "gamma" is the angle between the bond vectors.
    The distributions are stored if the "histogram" mode is selected.

    Parameters
    ----------
    universe : mda.Universe
        MDAnalysis universe. Only the bonds required for calculation of the
        ordering parameter shall be present. All other bonds shall be deleted
        from the universe before the analysis.
    r_min : FLOAT, optional
        Minimum distance between the bond vector centers to consider.
        The default is 0.. To exclude the bond itself, consider setting
        r_min to a small value, e. g. 1e-6
    r_max : FLOAT, optional
        Maximum distance between the bond vector centers to consider.
        The default is 0., which means choosing the cutoff based on the
        size of the simulation cell.
        Setting the value to -1, means considering all the bonds.
    mode : STRING, optional
        Whether an average value or a histogram shall be returned.
        The default is 'average'.
    same_molecule : BOOL, optional
        Whether the bonds from the same molecule (resid) shall be accounted
        for. The default is True.

    Returns
    -------
    Dictionary with "description" and "data" values. "data" contains the
    dictionaries for each timestep. For each timestep, the corresponding
    dictionary contains the key:value pairs of the parameters. The average
    value of the parameter s is calculated in the "average" mode, and
    the histograms normed by the total area in terms of cos^2(theta) and
    normed by the solid angle values are calculated in the "histogram" mode.
        
        {description: "Description containing the calculation parameters",
         data: {ts1: {"average_s": s, ...},
                ts2: {"average_s": s, ...},
                ....
               }
        }

    """
    # Prepare the description of the output:
    description = "Local orientational ordering parameter s"
    if mode == "histogram":
        description += (" and the distribution of " 
                     + "mutual orientation angle theta")
    description += (" r_min=" + str(r_min) + ", r_max=" + str(r_max))
    if same_molecule:
        description += ", same molecules taken into account"
    else:
        description += ", same molecules not taken into account"
    # Unwrap all the coordinates, so that all the bond lengths are
    # real. The closest images of the bonds will be found in the nested
    # function.
    unwrap = transformations.unwrap(u.atoms)
    u.trajectory.add_transformations(unwrap)
    
    data = {}
    
    for ts in u.trajectory:
        # Create data structure for the current timestep:
        values = {}
        #if mode == "average":
        cos_sq_sum = 0.
        i_s = 0
        if mode == "histogram":
            cos_sq_raw_hist = np.zeros(n_bins)
            _, bin_edges = np.histogram(cos_sq_raw_hist, bins = n_bins,
                                                       range = (0.,1.))
            # Calculate the values of cos_theta from cos_sq_theta:
            bin_edges_cosine = np.sqrt(bin_edges)
            # Calculate the values of theta corresponding to bin edges:
            solid_angle_normalization = np.diff(bin_edges_cosine)
            # Calculate the values of angle corresponding to bin edges:
            bin_edges_theta = np.arccos(bin_edges)
        
        if r_max == 0.:
            r_max = min(u.dimensions) / 2.
    
        # Calculate bond components
        # 1D arrays, one for each of the coordinates, provide more efficient
        # numpy calculations. Converting the data here, outside of the main
        # loop provided additional 15% speedup in the test runs.
        bx = (u.bonds.atom2.positions[:, 0] - u.bonds.atom1.positions[:, 0])
        by = (u.bonds.atom2.positions[:, 1] - u.bonds.atom1.positions[:, 1])
        bz = (u.bonds.atom2.positions[:, 2] - u.bonds.atom1.positions[:, 2])
    
        bond_components = [bx, by, bz]
    
        # Creating 1D arrays with bond midpoints
        rx = (u.bonds.atom1.positions[:, 0] 
              + u.bonds.atom2.positions[:, 0]) / 2.
        ry = (u.bonds.atom1.positions[:, 1]
              + u.bonds.atom2.positions[:, 1]) / 2.
        rz = (u.bonds.atom1.positions[:, 2]
              + u.bonds.atom2.positions[:, 2]) / 2.
    
        bond_midpoints = [rx, ry, rz]
    
        if not same_molecule:
            bond_resids = u.bonds.atom1.resids
        else:
            bond_resids = None
    
        for bond in u.bonds:
            # Determine the reference vector components and midpoint
            # from the bond coordinates
            ref_components = bond.atoms[1].position - bond.atoms[0].position
            ref_midpoint = (bond.atoms[0].position
                            + bond.atoms[1].position) / 2.
            # If needed, exclude bonds from the same molecule
            if not same_molecule:
                excluded_resids = bond.atoms[0].resid
            else:
                excluded_resids = None
            # Calculate ordering parameter value for the reference bond
            cos_sq_values = calculate_cos_sq_for_reference(
                bond_components, bond_midpoints, ref_components, ref_midpoint,
                u.dimensions, r_min = r_min, r_max = r_max,
                vector_attributes = bond_resids,
                excluded_attributes = excluded_resids)
            
            #if mode == "average":
            if np.shape(cos_sq_values)[0] > 0:
                cos_sq_sum += np.average(cos_sq_values)
                i_s += 1
            else:
                pass
                
            if mode == "histogram":
                cos_sq_hist_increment, _ = np.histogram(cos_sq_values,
                                           bins = n_bins, range = (0.,1.))
                cos_sq_raw_hist += cos_sq_hist_increment

        #if mode == "average":
        if i_s > 0:
            # Normalize the values. Normalization procedure ensures that
            # double consideration of each of the bonds doesn't affect
            # the result
            values["average_s"] = 1.5 * cos_sq_sum / i_s - 0.5
        if mode == "histogram":
            values["cos_sq_raw_histogram"] = cos_sq_raw_hist.tolist()
            norm = np.sum(cos_sq_raw_hist * np.diff(bin_edges))
            values["cos_sq_area_normalized_histogram"] = ( cos_sq_raw_hist
                                                            / norm).tolist()
            solid_angle_norm = np.sum(cos_sq_raw_hist * np.diff(bin_edges)
                                      / solid_angle_normalization )
            values["cos_sq_solid_angle_normalized_histogram"] = (
                cos_sq_raw_hist / solid_angle_normalization
                / solid_angle_norm ).tolist()
            values["bin_edges_cos_sq_theta"] = bin_edges.tolist()
            values["bin_edges_cos_theta"] = bin_edges_cosine.tolist()
            values["bin_edges_theta"] = bin_edges_theta.tolist()
        data[str(ts)] = values
    return { "description" : description, "data" : data }

def main():
    
    import argparse

    parser = argparse.ArgumentParser(
        description = """This utility calculates the angles between the bonds,
        if their midpoints are located within the range of [rmin, rmax].
        The local ordering parameter is then calculated as
        S = 3/2 ( <(cos(gamma))^2>) - 1/2)
        where "gamma" is the angle between the bond vectors.
        The distributions are stored if the --histogram flag is provided.
        The example applications are
        https://doi.org/10.1016/j.polymer.2020.122232
        and https://doi.org/10.1016/j.polymer.2022.124974""")

    parser.add_argument(
        'input', metavar = 'INPUT', action = "store", nargs = '+',
        help = """input file(s), the format will be guessed by MDAnalysis 
        based on file extension""")

    parser.add_argument(
        '--r_max', metavar = 'R_max', type = float, nargs = '?',
        default = 0., help = "outer cutoff radius")

    parser.add_argument(
        '--r_min', metavar = 'R_min', type = float, nargs = '?',
        default = 0., help = "inner cutoff radius")
    
    parser.add_argument(
        "--same-molecule", action = "store_true",
        help = "Take into account bonds from the same molecule")

    parser.add_argument(
        '--histogram', action = "store_true",
        help = "Store and optionally plot the distribution of the angles")

    parser.add_argument(
        '--n_bins', metavar = 'N_bins', type = int, nargs = '?',
        default = 150, help = "Number of bins of the distribution histogram")
    
    parser.add_argument('--plot', action = "store_true",
                        help = "Plot the distribution histogram")
    

    args = parser.parse_args()

    u = mda.Universe(*args.input)
    
    mode = "average"
    if args.histogram:
        mode = "histogram"
    
    result = calculate_orientation_order_parameter(u, r_min = args.r_min,
                                                   r_max = args.r_max,
                                                   mode = mode,
                                                   n_bins = args.n_bins,
                                                   same_molecule
                                                   = args.same_molecule,
                                                  )

    print(json.dumps(result, indent = 2))
    
    # Plot the histogram, if requested, with
    # the values summed across the timesteps
    if args.plot:
        import matplotlib.pyplot as plt
        frequencies, bincenters = averaged_frequencies_bin_centers(
          result, "cos_sq_area_normalized_histogram", "bin_edges_cos_sq_theta")
        plt.plot(bincenters, frequencies, label = "Histogram area normalized")
        frequencies, bincenters = averaged_frequencies_bin_centers(result,
           "cos_sq_solid_angle_normalized_histogram", "bin_edges_cos_sq_theta")
        plt.plot(bincenters, frequencies, label = "Solid angle normalized")
        
        plt.xlim(0, 1)
        plt.ylim(0)
        plt.yticks([0], fontsize = 18)
        plt.xlabel('cos_sq_\u03B8', fontsize = 18)
        plt.ylabel('P(\u03B8), a.u.', fontsize = 18)
        plt.legend(shadow = False, fontsize = 18)
        plt.show()
        
if __name__ == "__main__":
    main()