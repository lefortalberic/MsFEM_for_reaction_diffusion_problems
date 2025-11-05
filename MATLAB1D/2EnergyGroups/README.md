# Getting started
This folder corresponds to the simulation of the reaction-diffusion eigenvalue problem. It implements a MsFEM-type solution method for a vectorial equation in one dimension.

The main file to run the simulation is [main_diffusive_eigenproblem_ueps_vectoriel.m](main_diffusive_eigenproblem_ueps_vectoriel.m). The idea number (num_idee) corresponds to the type of basis function to be constructed in the MsFEM approximation space.
The method described in the article "Multiscale Finite
Element Method for reaction-diffusion eigenproblems and applications." is the Filter Method, i.e., method number 12.

The main_compare_* files allow plotting error comparisons.
