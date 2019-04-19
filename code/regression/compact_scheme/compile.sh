 #!/bin/bash

gfortran -O3 -fopenmp -ffree-line-length-512 calc_derivative.f90 -o calc_derivative.out
