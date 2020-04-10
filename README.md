# Quantum Dynamics in MATLAB

In this repository, we release the codes used for the paper:

https://iopscience.iop.org/article/10.1088/1361-6404/ab8360

This project uses MATLAB to simulate closed and open quantum systems. For the closed quantum system, we focus on the time-independent Schrödinger equation applied to the Ising model and cavity QED arrays. For open quantum systems, we introduce a fast algorithm to solve the Markovian master equation. Also, we address the Open Ising Model and a Two-Level-System coupled to a photon reservoir. Finally, we introduce the time-local master equation applied to the spin-boson model with pure dephasing. There are various functions used to optimize the presented dynamics.

# Functions

getSci.m: Generates the many-body Pauli matrices for S=1/2 and N particles

acav.m:  Generates the annihilation boson operator a for a system of N interacting cavities

sigmap.m:  Generates the atom operator sigma_+ a for a system of N interacting cavities

sortingEigenvalues.m:  Find the sorted eigenvalues and eigenmatrices of the Lindblad super-operator

QuantumSimulationCavityArray.m:  Solves the closed dynamics of cavity QED arrays for the Jaynes-Cumming-Hubbard and Rabi-Hubbard models

# Closed quantum systems

ClosedIsingModelTwoSpins.m:  Solves the Schrödinger dynamics of two interacting spins coupled to an external magnetic field

CavityQEDTransitionPhase.m:  Solves the Schrödinger dynamics and computes de Loschmidt echo for  a cavity QED array system

TransitionPhaseIsingModel.m:  Solves the Schrödinger dynamics and computes de Loschmidt echo for the transverse Ising model

# Open quantum systems

OpenIsingModelTwoSpins.m: Solves the quantum master equation for two interacting spins coupled to the same bosonic reservoir

TwoLevelSystemCoupledToLightBath.m: Solves the quantum master equation for a two-level system coupled to a photonic reservoir at finite temperature

NonMarkovianDynamicsPureDephasing.m: Solves the time-local quantum master equation (non-Markovian dynamics) for the pure-dephasing spin-boson model
