# Introduction

This repository contains all the code used to generate data for the article "Local equations describe unreasonably efficient stochastic algorithms in random K-SAT". The user can find programs to:
* Run the algorithms Focused Metropolis Search (FMS) and G-WalkSAT on random K-SAT instances.
* Perform the numerical integration of the Conditional Dynamic Approximation (CDA) with FMS and G-WalkSAT dynamic rules
* Perform the numerical integration of the Cavity Master Equation (CME) with FMS and G-WalkSAT dynamic rules
* Perform the numerical integration of the Dynamic Independent Neighbors Approximation (DINA) with FMS and G-WalkSAT dynamic rules
* Perform the numerical integration of the average version of CDA for FMS and G-WalkSAT dynamic rules

# Description of the different programs

## Algorithms

### FMS

File: wfacwsat.c
Parameters:
  The code can receive several parameters. For the article, we varied a subset of them. The script

