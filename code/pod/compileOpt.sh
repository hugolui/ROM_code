 #!/bin/bash

gfortran  -fautomatic -fmax-stack-var-size=8192  -fdefault-real-8  -fdefault-double-8 -ffree-line-length-none -O2  -I /usr/local/include -L /usr/local/lib -lcgns modfile.f90 cgns_routines.f90 compute_area.f90 pod.f90 -lcgns -o pod

