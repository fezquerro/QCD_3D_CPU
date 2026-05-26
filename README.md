
# Monte Carlo SU(2) Lattice Gauge Theory Parallel Simulator (3+1)


**A high-performance lattice simulation code for parallel Monte Carlo simulation for SU(2) Yang-Mills theory in 3+1 dimensions in C.**



## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Running](#running)

## Overview

This project implements a Monte Carlo simulation of pure SU(2) lattice gauge theory in four-dimensional Euclidean spacetime using the Wilson plaquette action. The average plaquette and internal energy are computed. The parallelization is applied on the temporal side of the lattice which allows an odd number of sites on the other direction of the lattice. This parallelization is applied using OpenMP. Two different boundary conditions are implemented for the third spatial dimensions, periodic and perfect colour conductor boundaty conditions. 

This code was primarily developed and used to generate the results presented in the PhD thesis  **[Properties of the quantum vacuum in non-abelian gauge theories](https://arxiv.org/abs/2605.18220)**. See Chapter 4 for an explanation of the two different boundary conditions and parallelization. The results of this code are presented in Chapter 6 and Appendix C.

## Features

- **Efficient SU(2) updates**: Heat-bath algorithm + over-relaxation
- **Flexible lattice sizes**: Support for different $N_t \times N_x \times N_y \times N_z$ lattices size and coupling constants $\beta$
- **Different boundary conditions**: Periodic and perfect colour conductor boundary conditions are implemented.
- **Parallelization**: OpenMP support for multi-core parallelization
- **Observables**:
  - Plaquette (average action)
  - Energy (first component of the energy momemtum tensor)
- **Thermalization and autocorrelation analysis**
- **Checkpointing** for long runs
- **Configurable random number generator** (with reproducible seeds)

## Running

For compiling in Linux 
```bash
gcc -o ejecutable main.c -lm -fopenmp
```
For compiling in Windows IDEs like Codeblocks add
``` bash
-fopenmp
```
in "other compiler options" in the compiler settings.

