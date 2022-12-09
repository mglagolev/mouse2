#!/usr/bin/env python3

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


import numpy as np
import sys
import MDAnalysis as mda
from MDAnalysis import transformations

def normalize_vectors(vectors):
    """
    Normalize an array of vectors

    """
    shape = vectors.shape
    norms = np.linalg.norm(vectors, axis = 1).reshape((shape[0], 1))
    return vectors / norms

def normal_vector(dir_vectors):
    """
    Calculate the normal vector, as the vector corresponding
    to the largest eigenvalue.

    """
    # compute the gyration tensor for the vectors connecting COMs
    tensor = (np.einsum('im,in->mn', dir_vectors, dir_vectors) 
                   / float(len(dir_vectors)))
    # find eigenvalues and eigenvectors 
    eigen_vals, eigen_vecs = np.linalg.eig(tensor)
    # the maximal eigenvalue
    max_col = list(eigen_vals).index(max(eigen_vals))
    # the corresponding eigenvector
    normal_vector = eigen_vecs[:,max_col]
    return normal_vector

def pk_v(vectors, reference_2, reference_3):
    """
    Calculate order parameters Pk and V,
    as described in https://doi.org/10.1063/5.0005854

    """
    c = np.cross(reference_2, reference_3)
    # The following operations are performed for arrays:
    sc = np.dot(vectors, reference_2)
    # Determine the angle between the vectors in the range [0,pi]
    gamma = np.arccos( sc / np.linalg.norm(vectors, axis = 1)
                          / np.linalg.norm(reference_2))
    sc_r = np.outer(sc, reference_2)
    rot_ak = vectors - sc_r
    phi = np.arccos( np.dot(rot_ak, c)
                    / np.linalg.norm(rot_ak, axis = 1)
                    / np.linalg.norm(c))
    pk = np.sin(gamma) * np.sin(gamma) * np.cos( 2. * phi)
    v = np.sin(2. * gamma) * np.cos(phi)
    average_pk = np.average(pk)
    average_v = np.average(v)
    return average_pk, average_v


def lamellar_ordering_parameters(u: mda.Universe, type_A, type_B):
    # Unwrap the atom coordinates
    unwrap = transformations.unwrap(u.atoms)
    u.trajectory.add_transformations(unwrap)
    for ts in u.trajectory:
        # Create data structures for the values for individual residues
        # Optimize performance by creating lists and then converting to NumPy
        A_start_list, A_end_list, A_com_list = [], [], []
        B_start_list, B_end_list, B_com_list = [], [], []
        # Create the list of unique residue ids:
        resids = set(u.atoms.resids)
        # For resid in residue id list:
        for resid in resids:
            # select atoms (resid, type1), select atoms (resid, type2)
            atoms_A = u.select_atoms("resid " + str(resid)
                                    + " and type " + str(type_A))
            atoms_B = u.select_atoms("resid " + str(resid)
                                    + " and type " + str(type_B))
            # The AtomGroups in MDAnalysis are ordered, so we can take
            # the 1st and the last atom of the groups:
            A_start_list.append(atoms_A[0].position)
            A_end_list.append(atoms_A[-1].position)
            A_com_list.append(atoms_A.center_of_mass())
            # Same for type B particles:
            B_start_list.append(atoms_B[0].position)
            B_end_list.append(atoms_B[-1].position)
            B_com_list.append(atoms_B.center_of_mass())
        # Convert lists to NumPy arrays:
        A_start = np.asarray(A_start_list)
        A_end = np.asarray(A_end_list)
        A_com = np.asarray(A_com_list)
        # Same for type B particles:
        B_start = np.asarray(B_start_list)
        B_end = np.asarray(B_end_list)
        B_com = np.asarray(B_com_list)
        # Calculate the end-to-end vectors for both blocks and the vectors
        # between the centers of mass of the blocks
        block_A_vectors_normed = normalize_vectors(A_end - A_start)
        block_B_vectors_normed = normalize_vectors(B_end - B_start)
        com_vectors_normed = normalize_vectors(B_com - A_com)
        lam_norm = normal_vector(com_vectors_normed)
        # the normal to the lamella
        sys.stderr.write("Eigen vector "+str(lam_norm)+"\n")
        sk_B = 1.5 * np.square(np.dot(block_B_vectors_normed, lam_norm)) - 0.5
        sk_B_average = np.average(sk_B)
        sk_B_histogram = np.histogram(sk_B)
        sys.stderr.write("Sk_B parameter = "+str(sk_B_average)+"\n")
        sys.stderr.write("Sk_B histogram = "+str(sk_B_histogram)+"\n")
        # calculate the directors for both blocks using the same approach
        block_A_director = normal_vector(block_A_vectors_normed)
        block_B_director = normal_vector(block_B_vectors_normed)
        sys.stderr.write("Eigen vector for B director "
                         + str(block_B_director)+"\n")
        h_B = np.cross(block_B_director, lam_norm)
        sys.stderr.write("h_B vector "+str(h_B)+"\n")
        if np.linalg.norm(h_B) == 0.:
            sys.stderr.write("Vector h_B = 0"+"\n")
        else:
            pk_B, v_B = pk_v(block_B_vectors_normed, lam_norm, h_B)
            tan_2theta = v_B / (sk_B_average - 0.5 * pk_B)
            theta = np.arctan (tan_2theta) / 2.0
            sys.stderr.write("P_k_B = " + str(pk_B) + "; V_B = "
                             + str(v_B) + "\n")
            sys.stderr.write("tan_2theta = " + str(tan_2theta) 
                             + "; theta = " + str(theta) + "\n")
        
if __name__ == "__main__":
    
    import argparse

    parser = argparse.ArgumentParser(
        description = 'Calculate the lamellar ordering parameter')

    parser.add_argument(
        'input', metavar = 'INPUT', action = "store", nargs = '+',
        help = "input files")

    parser.add_argument(
        '--block-types', metavar = 'TYPES', type = str, nargs = 2,
        default = ["1", "2"],
        help = "bead types for the blocks (provide 0 or 2 arguments)")
    

    args = parser.parse_args()

    u = mda.Universe(*args.input)
    
    lamellar_ordering_parameters(u, args.block_types[0],
                                             args.block_types[1])