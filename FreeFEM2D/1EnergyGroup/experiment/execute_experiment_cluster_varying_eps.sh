#! /bin/bash

# To be executed from the directory where main_LIN_MPI.edp and main_CR_MPI.edp are located

# Consecutive executions of FreeFem++ code to perform experiments with different parameter values
# The parameter values are changed continuously based on the values given below
# Check that all values that are not controlled here, are as desired in the basic parameter.txt file in experiment_parameters
# The execution commands for FreeFem++ are adapted to use on the cluster of CERMICS


# Number of processors to be used
NUMBER_OF_PROC=1

# Parameter values to be used in the tests (all will be combined)
# eg TOTEST_LARGE_N="8 16 32" to test for three different (coarse) mesh sizes
TOTEST_L="1."
TOTEST_LARGE_N="1800" 
TOTEST_SMALL_N="3"
#TOTEST_EPS="0.1667 0.148 0.1250 0.117 0.1 0.093 0.0887 0.0836 0.08 0.0779 0.0738 0.0697 0.0656 0.0625 0.0575 0.0543 0.05 0.046 0.038 0.033"
#TOTEST_EPS="0.1667 0.125 0.1 0.0887 0.0779 0.0675 0.0575 0.04545 0.04032 0.033333 0.02873 0.026316"
TOTEST_EPS="0.6668 0.61 0.53 0.47 0.4 0.36 0.3116 0.27 0.23 0.21 0.1818 0.16128 0.133332 0.11492 0.105264 0.092 0.073 0.054"
#TOTEST_EPS="0.6668 0.61 0.53 0.47 0.4 0.36 0.3116 0.27 0.23 0.21 0.1818"
#TOTEST_EPS="0.054"
TOTEST_COEFF="2"
TOTEST_VP="1"
#TOTEST_OS="1 1.199988 1.466652 1.733316 2 2.2666 2.5333 2.7999 3.0666 3.3333 4.133292 4.933284"
TOTEST_OS="2"


# LOOP OVER ALL PARAMETER VALUES AND FreeFem++ EXECUTION
for TEST_COEFF in $TOTEST_COEFF; do sed -i "s/NumCoeffDiff=.*/NumCoeffDiff= $TEST_COEFF/" "experiment/parameters.txt"
for TEST_VP in $TOTEST_VP; do sed -i "s/NumberEigenValue=.*/NumberEigenValue= $TEST_VP/" "experiment/parameters.txt"
for TEST_L in $TOTEST_L; do sed -i "s/L=.*/L= $TEST_L/" "experiment/parameters.txt" 
for TEST_EPS in $TOTEST_EPS; do sed -i "s/eps=.*/eps= $TEST_EPS/" "experiment/parameters.txt" 
for TEST_LARGE_N in $TOTEST_LARGE_N; do sed -i "s/N=.*/N= $TEST_LARGE_N/" "experiment/parameters.txt"
    # The above loops contain all parameters related to the reference solution

    # Reference solution
    cp experiment/parameters.txt parameters.txt


    for TEST_SMALL_N in $TOTEST_SMALL_N; do sed -i "s/n=.*/n= $TEST_SMALL_N/" "experiment/parameters.txt" 
    # The above loops contain all parameters that require the computation of a new basis

        cp experiment/parameters.txt parameters.txt
        
        #mpirun -np $NUMBER_OF_PROC FreeFem++-mpi main_react_diff_triche.edp -v 0
        mpirun -np $NUMBER_OF_PROC FreeFem++-mpi main_react_diff_P1.edp -v 0
        
        for TEST_OS in $TOTEST_OS; do sed -i "s/osCoef=.*/osCoef= $TEST_OS/" "experiment/parameters.txt"
            
            cp experiment/parameters.txt parameters.txt

            #mpirun -np $NUMBER_OF_PROC FreeFem++-mpi main_react_diff_MsFEM_OS_square.edp -v 0
            mpirun -np $NUMBER_OF_PROC FreeFem++-mpi main_react_diff_MsFEM_OS_square_filtre.edp -v 0
            #mpirun -np $NUMBER_OF_PROC FreeFem++-mpi main_react_diff_MsFEM_OS_square_filtre_analysis.edp -v 0
            #mpirun -np $NUMBER_OF_PROC FreeFem++-mpi main_erreur_phi_psi_cellule_unite_filtre_recombine.edp -v 0

        done

    done # end of loops over numerical parameters
done done done done done # end of loops over reference solution/physical parameters + fine mesh
