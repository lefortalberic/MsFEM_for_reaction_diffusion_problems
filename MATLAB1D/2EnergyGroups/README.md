# Getting started
This folder corresponds to the simulation of the reaction-diffusion eigenvalue problem. It implements a MsFEM-type solution method for a vectorial equation in one dimension.

The main file to run the simulation is [main_diffusive_eigenproblem_ueps_vectoriel.m](main_diffusive_eigenproblem_ueps_vectoriel.m). The idea number (num_idee) corresponds to the type of basis function to be constructed in the MsFEM approximation space.

The approach presented in the article by **Claude Le Bris, Albéric Lefort, and Frédéric Legoll**,  
*Multiscale Finite Element Method for Reaction–Diffusion Eigenproblems and Applications*,  
corresponds to the **Filter Method** (method no. 12).

The main_compare_* files allow plotting error comparisons.
