# Introduction

This repository contains all the code used to generate data for the article "Local equations describe unreasonably efficient stochastic algorithms in random K-SAT" ([arxiv](http://arxiv.org/abs/2504.06757)). The user can find programs to:
* Run the algorithms Focused Metropolis Search (FMS) and G-WalkSAT on random K-SAT instances.
* Perform the numerical integration of the Conditional Dynamic Approximation (CDA) with FMS and G-WalkSAT dynamic rules
* Perform the numerical integration of the Cavity Master Equation (CME) with FMS and G-WalkSAT dynamic rules
* Perform the numerical integration of the Dynamic Independent Neighbors Approximation (DINA) with FMS and G-WalkSAT dynamic rules
* Perform the numerical integration of the average version of CDA for FMS and G-WalkSAT dynamic rules
* Run the CDA-guided decimation with FMS dynamic rules

It also contains most of the data necessary to reobtain the figures in the paper, and the scripts that produce them

# Description of the different programs

## Algorithms

### FMS

File: wfacwsat.c
This implements the Focused Metropolis Search (FMS) algorithm for solving the **K-SAT problem** on random graphs [[1]](#1). This is an implementation by Sakari Seitz, made as a modification of the public walksat package, by Henry Kautz.
Auxiliary files: Graph_to_CNF_input.cpp, Scripts/run_FMS_KSAT_PD_inner.sh
Parameters:
  The code can receive several parameters. For the article, we varied a subset of them using the script "Scripts/run_FMS_KSAT_PD_inner.sh". The script receives:

  * N -> number of variables in the formula
  * alpha -> ratio between the number of clauses and the number of variables
  * K -> variables per clause
  * nsamples -> number of different formulas (instances) to generate and attempt to solve
  * hist -> number of tries per instance
  * seed_fms -> seed for the random number generator
  * iters_save -> the program prints the energy every 'iters_save' iterations
  * tl -> maximum time to run, in Monte Carlo sweeps
  * eta -> algorithmic parameter. In the algorithm, eta is real number between zero and one. This script receives an integer to be divided by 100 internally in the code to get the right value of eta.
  * path -> path to the output files
  * print_every -> the program saves statistics every 'print_every' iterations to a file.

Language: C, C++
Requires: GSL (GNU Scientific Library)


### G-WalkSAT

File: walksat_greedy.JL
Auxiliary file: Scripts/run_walksat-greedy_KSAT_PD_1.sh (for an example of how to run the code)
Parameters:
  The julia script 'walksat_greedy.JL' receives:

  * N -> number of variables in the formula
  * M -> number of clauses in the formula
  * q -> algorithmic parameter of G-WalkSAT
  * K -> variables per clause
  * Exp -> number of different formulas (instances) to generate and attempt to solve
  * Hist -> number of tries per instance
  * time -> maximum time to run, in Monte Carlo sweeps
  * ruta -> path to the output files
  * save_every -> the program saves statistics every 'print_every' iterations to a file.


Language: julia
Requires: SimpleHypergraphs.jl (is a julia package)


## CDA single instance

### CDA for FMS

File: CDA_FMS.cpp
Auxiliary file: Scripts/run_CDA_cpp.sh (for an example of how to run the code)
Parameters:
  The cpp file 'CDA_FMS.cpp' receives the command line parameters:

  * N -> number of variables in the formula
  * M -> number of clauses in the formula
  * K -> variables per clause
  * seed_r -> seed for the generation of the boolean formula
  * eta -> FMS's algorithmic parameter
  * tl -> time limit in Monte Carlo sweeps
  * tol -> tolerance for the numerical integration of the differential equations
  * nthr -> number of threads to use during execution.

Language: C++
Requires: GSL (GNU Scientific Library), OpenMP

### CDA for G-WalkSAT

File: CDA_WalkSAT_av_rates.cpp
Auxiliary file: Scripts/run_CDA_WalkSAT_av_rates.sh (for an example of how to run the code)
Parameters:
  The cpp file 'CDA_WalkSAT_av_rates.cpp' receives the command line parameters:

  * N -> number of variables in the formula
  * M -> number of clauses in the formula
  * K -> variables per clause
  * seed_r -> seed for the generation of the boolean formula
  * q -> G-WalkSAT's algorithmic parameter
  * tl -> time limit in Monte Carlo sweeps
  * tol -> tolerance for the numerical integration of the differential equations
  * nthr -> number of threads to use during execution.

Language: C++
Requires: GSL (GNU Scientific Library), OpenMP


## CME single instance

### CME for FMS

File: CME_FMS.cpp
Auxiliary file: Scripts/run_CME_cpp.sh (for an example of how to run the code)
Parameters:
  The cpp file 'CME_FMS.cpp' receives the command line parameters:

  * N -> number of variables in the formula
  * M -> number of clauses in the formula
  * K -> variables per clause
  * seed_r -> seed for the generation of the boolean formula
  * eta -> FMS's algorithmic parameter
  * tl -> time limit in Monte Carlo sweeps
  * tol -> tolerance for the numerical integration of the differential equations
  * nthr -> number of threads to use during execution.

Language: C++
Requires: GSL (GNU Scientific Library), OpenMP


### CME for G-WalkSAT

File: CME_WalkSAT_av_rates.cpp
Auxiliary file: Scripts/run_CME_WalkSAT_av_rates.sh (for an example of how to run the code)
Parameters:
  The cpp file 'CME_WalkSAT_av_rates.cpp' receives the command line parameters:

  * N -> number of variables in the formula
  * M -> number of clauses in the formula
  * K -> variables per clause
  * seed_r -> seed for the generation of the boolean formula
  * q -> G-WalkSAT's algorithmic parameter
  * tl -> time limit in Monte Carlo sweeps
  * tol -> tolerance for the numerical integration of the differential equations
  * nthr -> number of threads to use during execution.

Language: C++
Requires: GSL (GNU Scientific Library), OpenMP


## DINA

### DINA for FMS

File: Weigt_KSAT_cluster.py
Auxiliary file: Scripts/scripts_Weigt_KSAT_FMS_tasks_1.sh (for an example of how to run the code)
Parameters:
  The python script 'Weigt_KSAT_cluster.py' receives the command line parameters:

  * eta -> FMS's algorithmic parameter
  * alpha -> ratio between the number of clauses and the number of variables
  * tl -> time limit in Monte Carlo sweeps

Language: Python
Requires: numpy, scipy


### DINA for G-WalkSAT

File: DINA_G-WalkSAT.JL
Parameters:
  The julia script 'DINA_G-WalkSAT.JL' receives the command line parameters:

  * q -> G-WalkSAT's algorithmic parameter

The values of the rest of the parameters can be adjusted inside the script

Language: Julia
Requires: DifferentialEquations.jl


## Average case CDA

### Average case CDA for FMS

File: CDA1av_FMS_lpln_pop_dyn.cpp
Auxiliary file: Scripts/run_CDA1av_lpln_pop_dyn_K_3.sh (for an example of how to run the code)
Parameters:
  The cpp file 'CDA1av_FMS_lpln_pop_dyn.cpp' receives the command line parameters:

  * pop_size -> size of the population of probabilities (see Supporting Information Appendix for the article)
  * alpha -> ratio between the number of clauses and the number of variables
  * K -> variables per clause
  * seed_r -> seed for the generation of the boolean formula
  * eta -> FMS's algorithmic parameter
  * tl -> time limit in Monte Carlo sweeps
  * tol -> tolerance for the numerical integration of the differential equations
  * nthr -> number of threads to use during execution.
  * eps_c -> initially, the program selects the maximum connectivity in the population by finding
  the first c_max such that Poisson(cmax, <c>) < eps_c 

Language: C++
Requires: GSL (GNU Scientific Library), OpenMP


### Average case CDA for G-WalkSAT

File: CDA1av_WalkSAT_lpln_pop_dyn.cpp
Auxiliary file: Scripts/run_CDA1av_lpln_WalkSAT_pop_dyn_K_3_1.sh (for an example of how to run the code)
Parameters:
  The cpp file 'CDA1av_WalkSAT_lpln_pop_dyn.cpp' receives the command line parameters:

  * pop_size -> size of the population of probabilities (see Supporting Information Appendix for the article)
  * alpha -> ratio between the number of clauses and the number of variables
  * K -> variables per clause
  * seed_r -> seed for the generation of the boolean formula
  * q -> G-WalkSAT's algorithmic parameter
  * tl -> time limit in Monte Carlo sweeps
  * tol -> tolerance for the numerical integration of the differential equations
  * nthr -> number of threads to use during execution.
  * eps_c -> initially, the program selects the maximum connectivity in the population by finding
  the first c_max such that Poisson(cmax, <c>) < eps_c 

Language: C++
Requires: GSL (GNU Scientific Library), OpenMP


## CDA-guided decimation with FMS's rates

File: CDA_decimation_FMS.cpp
Auxiliary file: Scripts/run_CDA_decimation_FMS_inside.sh (for an example of how to run the code)
Parameters:
  The cpp file 'CDA_decimation_FMS.cpp' receives the command line parameters:

  * N -> number of variables in the formula
  * M -> number of clauses in the formula
  * K -> variables per clause
  * seed_r -> seed for the generation of the boolean formula
  * eta -> FMS's algorithmic parameter
  * steps_dec -> number of integrator steps between consecutive decimations.
  * tol -> tolerance for the numerical integration of the differential equations
  * nthr -> number of threads to use during execution.

Language: C++
Requires: GSL (GNU Scientific Library), OpenMP

# Data

The folder Data/ contains most of the data necessary to reobtain the figures in our paper ([arxiv](http://arxiv.org/abs/2504.06757)). Inside, the user can find a folder for each one of the figures. Some of the data could not be uploaded to GitHub because of space issues (in total, it would be more than 30GB).

# References
<a id="1">[1]</a> 
Sakari Seitz, Mikko Alava and Pekka Orponen. 
Focused local search for random 3-satisfiability. 
J. Stat. Mech. Teory Exp., P06006 (2005).