Code in C for computing in parallel the average value of the plaquette and the energy in a 3+1 SU(2) lattice.

The parallelization is implemented using OpenMP https://www.openmp.org/. The temporal side of the lattice is the one parallelized which allows an odd number of sites on the other directions of the lattice.

For compiling in Linux 
```bash
gcc -o ejecutable main.c -lm -fopenmp
```
For compiling in Windows IDEs like Codeblocks add
``` bash
-fopenmp
```
in "other compiler options" in the compiler settings.
