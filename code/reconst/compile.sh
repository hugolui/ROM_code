 #!/bin/bash

gfortran  -fautomatic -fmax-stack-var-size=8192  -fdefault-real-8  -fdefault-double-8 -ffree-line-length-none -O2 -I /usr/local/include -L /usr/local/lib -lcgns modfile.f90 reconst.f90 -lcgns -o reconst.out

