# MsFEM
In this project we implement Multiscale Finite Element Method (MsFEM) type approaches for
reaction–diffusion equations with oscillatory coefficients. These approaches are based on a Galerkin-type formulation on a coarse mesh.
The basis functions of the approximation space are solutions to local problems, and are thus well adapted to the multiscale problem under consideration.
These functions are constructed locally on each element of the coarse mesh using a mesh sufficiently fine to capture
all the details of the microstructure.

The implementation of MsFEM for reaction-diffusion eigenproblems are in [MATLAB](https://fr.mathworks.com/) for 1D simulations and in [FreeFEM](https://freefem.org/) for 2D simulations.
This project is so divided into the two different parts detailed below.

Each part contains a ReadMe file detailing the operation of the respective part.

The main approach is presented in the article by **Claude Le Bris, Albéric Lefort, and Frédéric Legoll**,  
*Multiscale Finite Element Method for Reaction–Diffusion Eigenproblems and Applications*.

## 1D Simulations
The first part is coded in MATLAB for all 1D simulations, providing a more accessible approach and manipulation.
This part is further divided into two sections. The first corresponds to the simulation of the scalar eigenvalue problem modeling a single energy group (in the [MATLAB1D/1EnergyGroup](MATLAB1D/1EnergyGroup) folder).
The second corresponds to the simulation of the vector eigenvalue problem modeling two energy groups (in the [MATLAB1D/2EnergyGroups](MATLAB1D/2EnergyGroups) folder).

More information about the 1D code can be found in the [wiki](https://github.com/lefortalberic/MsFEM_for_reaction_diffusion_problems/wiki/MATLAB-code-explanation) page.

## 2D Simulations
The second part corresponds to 2D simulations and is coded in FreeFEM. This part is also divided into two sections. The first corresponds to the simulation of the scalar eigenvalue problem modeling a single energy group (in the [FreeFEM2D/1EnergyGroup](FreeFEM2D/1EnergyGroup) folder).
The second corresponds to the simulation of the vector eigenvalue problem modeling two energy groups (in the [FreeFEM2D/2EnergyGroups](FreeFEM2D/2EnergyGroups) folder).

These simulations require the SLEPc package (Scalable Library for Eigenvalue Problem Computations) to solve eigenvalue problems efficiently.
Make sure that both [PETSc](https://petsc.org/) and [SLEPc](https://slepc.upv.es/) are correctly installed and configured with your FreeFEM environment.

More information about the 2D code can be found in the [wiki](https://github.com/lefortalberic/MsFEM_for_reaction_diffusion_problems/wiki/FreeFem-code-explanation) page.

# Contributing
Feel free to submit issues or pull requests if you would like to contribute to this project.
