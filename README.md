# mouse2
Scripts for molecular ordering analysis - new version based on NumPy and MDAnalysis.

This repository contains the utilities to quantitatively assess the results of molecular dynamics simulations by calculating numeric ordering parameters.
The focus is the chirality of the systems and the local ordering arising from it.

The NumPy and MDAnalysis libraries need to be installed to use the scripts, and the
NetworkX library is additionaly required to use the aggregates.py
Matplotlib and SciPy are required to use the plotting and fitting options in some of the scripts.

Quick reference:

usage: backbone_bond_autocorrelations.py [-h] [--k_max [k_max]] [--selection [QUERY]] [--different-molecules] [--plot] [--fit] [--p_guess [NUMBER]] INPUT [INPUT ...]


Calculate the autocorrelation function of the polymer bonds.
The formula is presented in https://doi.org/10.1134/S0965545X10070102
Application to all-atom simulations: https://doi.org/10.3390/polym11122056


positional arguments:
  INPUT                input file(s), the format will be guessed by MDAnalysis based on file extension

options:
  -h, --help            show this help message and exit
  --k_max [k_max]       maximum distance between the bonds along the backbone
  --selection [QUERY]   Consider only selected atoms, use MDAnalysis selection language
  --different-molecules
                        Calculate correlations based on particle index number, even if the bonds belong to different molecules
  --plot                Plot the averaged results
  --fit                 Fit the averaged results with a modulated exponential function
  --p_guess [NUMBER]    Initial guess for the number of monomer units per turn
                        


bond_orientational_ordering.py [-h] [--r_max [R_max]] [--r_min [R_min]] [--same-molecule] [--histogram] [--n_bins [N_bins]] [--plot] INPUT [INPUT ...]

This utility calculates the angles between the bonds, if their midpoints are located within the range of [rmin, rmax].
The local ordering parameter is then calculated as S = 3/2 ( <(cos(gamma))^2>) -1/2)
where "gamma" is the angle between the bond vectors. The distributions are stored if the --histogram flag is provided.
The example applications are https://doi.org/10.1016/j.polymer.2020.122232
and https://doi.org/10.1016/j.polymer.2022.124974



positional arguments:
  INPUT              input file(s), the format will be guessed by MDAnalysis based on file extension

options:
  -h, --help         show this help message and exit
  --r_max [R_max]    outer cutoff radius
  --r_min [R_min]    inner cutoff radius
  --same-molecule    Take into account bonds from the same molecule
  --histogram        Store and optionally plot the distribution of the angles
  --n_bins [N_bins]  Number of bins of the distribution histogram
  --plot             Plot the distribution histogram


lamellar_orientational_ordering.py [-h] [--block-types TYPES TYPES] [--A] [--B] [--verbose] INPUT [INPUT ...]

Calculate the molecular ordering parameters for lamellae containing tilted copolymer blocks, as described in the paper by 
M. A. Osipov, M. V. Gorkunov, A. V. Berezkin, A. A. Antonov and Y. V. Kudryavtsev
"Molecular theory of the tilting transition and computer simulations of the tilted lamellar phase of rod???coil diblock copolymers"
https://doi.org/10.1063/5.0005854.
A use case is also presented in https://doi.org/10.1039/D1SM00759A


positional arguments:
  INPUT                 input file(s), the format will be guessed by MDAnalysis based on file extension

options:
  -h, --help            show this help message and exit
  --block-types TYPES TYPES
                        bead types for the blocks A dnd B (provide 2 arguments, without the option default values 1 and 2 are used)
  --A                   Calculate the values for block A
  --B                   Calculate the values for block B
  --verbose             Store the values for individual molecules
  


backbone_twist.py [-h] [--selection [QUERY]] [--k VECTOR_LENGTHS [VECTOR_LENGTHS ...]] [--different-molecules] [--plot] INPUT [INPUT ...]

Calculate the list of dihedral angles, formed by the following vectors:
(r_i, r_{i+k}), (r_{i+k}, r_{i+2*k}), (r_{i+2*k}, r_{i+3*k})
where i is the index of a monomer unit.
The example of the analysis is provided in the Supplementary Information for
https://doi.org/10.1016/j.polymer.2022.124974


positional arguments:
  INPUT                 input file(s), the format will be guessed by MDAnalysis based on file extension

options:
  -h, --help            show this help message and exit
  --selection [QUERY]   Consider only selected atoms, use MDAnalysis selection language
  --k VECTOR_LENGTHS [VECTOR_LENGTHS ...]
                        List of vector lengths along the backbone
  --different-molecules
                        Consider the angles spanning different molecules
  --plot                Plot the results
  
  

aggregates.py [-h] [--r_neigh [R_neigh]] [--selection [QUERY]] INPUT [INPUT ...]

This utility returns a data structure containing list of aggregates for all of the timesteps in the MDAnalysis universe.
Each aggregate is determined as a complete graph of neighbors.
The atoms are considered neighbors if the distance between their centers does not exceed r_neigh.
Each aggregate is represented as a list of MDAnalysis atom indices.

positional arguments:
  INPUT                input file(s), the format will be guessed by MDAnalysis based on file extension

options:
  -h, --help           show this help message and exit
  --r_neigh [R_neigh]  neighbor cutoff
  --selection [QUERY]  Consider only selected atoms, use MDAnalysis selection language



The algorithms were used in the following publications:

Abramova A. A., Glagolev M. K., Vasilevskaya V. V. Structured globules with twisted arrangement of helical blocks: Computer simulation // Polymer. ??? 2022. ??? P. 124974.

Glagolev M. K., Glagoleva A. A., Vasilevskaya V. V. Microphase separation in helix-coil block copolymer melts: computer simulation // Soft Matter. ??? 2021. ??? Vol. 17, no. 36. ??? P. 8331???8342.

Glagolev M. K., Vasilevskaya V. V. Coarse-grained simulation of molecular ordering in polylactic blends under uniaxial strain // Polymer. ??? 2020. ??? Vol. 190. ??? P. 122232.

Glagolev M. K., Vasilevskaya V. V. Liquid-crystalline ordering of filaments formed by bidisperse amphiphilic macromolecules // Polymer Science - Series C. ??? 2018. ??? Vol. 60, no. 1. ??? P. 39???47.

Glagolev M. K., Vasilevskaya V. V., Khokhlov A. R. Domains in mixtures of amphiphilic macromolecules with different stiffness of backbone // Polymer. ??? 2017. ??? Vol. 125. ??? P. 234???240.

Glagolev M. K., Vasilevskaya V. V., Khokhlov A. R. Induced liquid-crystalline ordering in solutions of stiff and flexible amphiphilic macromolecules: Effect of mixture composition // Journal of Chemical Physics. ??? 2016. ??? Vol. 145, no. 4. ??? P. 044904.

Glagolev M. K., Vasilevskaya V. V., Khokhlov A. R. Effect of induced self-organization in mixtures of amphiphilic macromolecules with different stiffness // Macromolecules. ??? 2015. ??? Vol. 48, no. 11. ??? P. 3767???3774.

Glagolev M. K., Vasilevskaya V. V., Khokhlov A. R. Self-organization of amphiphilic macromolecules with local helix structure in concentrated solutions // Journal of Chemical Physics. ??? 2012. ??? Vol. 137, no. 8.

Glagolev M. K., Vasilevskaya V. V., Khokhlov A. R. Formation of fibrillar aggregates in concentrated solutions of rigid-chain amphiphilic macromolecules with fixed torsion and bend angles // Polymer Science, Series A. ??? 2011. ??? Vol. 53, no. 8. ??? P. 733???743.

Glagolev M. K., Vasilevskaya V. V., Khokhlov A. R. Compactization of rigid-chain amphiphilic macromolecules with local helical structure // Polymer Science, Series A. ??? 2010. ??? Vol. 52, no. 7. ??? P. 761???774.
