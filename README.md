# MsFEM
Implementation of MsFEM for reaction-diffusion eigenproblems in [FreeFEM](https://freefem.org/) and [MATLAB](https://fr.mathworks.com/).
This project is divided into two different parts detailed below.

Each part contains a ReadMe file detailing the operation of the respective part.

## 1D Simulations
The first part is coded in MATLAB for all 1D simulations, providing a more accessible approach and manipulation.
This part is further divided into two sections. The first corresponds to the simulation of the scalar eigenvalue problem modeling a single energy group (in the [MATLAB1D/1EnergyGroup](MATLAB1D/1EnergyGroup) folder).
The second corresponds to the simulation of the vector eigenvalue problem modeling two energy groups (in the [MATLAB1D/2EnergyGroups](MATLAB1D/2EnergyGroups) folder).

## 2D Simulations
The second part corresponds to 2D simulations and is coded in FreeFEM. This part is also divided into two sections. The first corresponds to the simulation of the scalar eigenvalue problem modeling a single energy group (in the [FreeFEM2D/1EnergyGroup](FreeFEM2D/1EnergyGroup) folder).
The second corresponds to the simulation of the vector eigenvalue problem modeling two energy groups (in the [FreeFEM2D/2EnergyGroups](FreeFEM2D/2EnergyGroups) folder).

These simulations require the SLEPc package (Scalable Library for Eigenvalue Problem Computations) to solve eigenvalue problems efficiently.
Make sure that both [PETSc](https://petsc.org/) and [SLEPc](https://slepc.upv.es/) are correctly installed and configured with your FreeFEM environment.

# Contributing
Feel free to submit issues or pull requests if you would like to contribute to this project.
