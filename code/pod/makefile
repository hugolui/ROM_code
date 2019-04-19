#=======================================================================
# User definitions
#

NAME = pod              

F90FILES = modfile.f90 pod.f90 cgns_routines.f90 compute_area.f90 CorrMatrix.f90 

# Ifort stuff
# ARCH = -heap-arrays -D_REENTRANT -r8 -convert big_endian  -L$(MKLPATH) -I$(MKLINCLUDE) -B$(PWD)/../libs 
# MKLPATH = /usr/local/intel2013/mkl/lib/intel64
# MKLINCLUDE = /usr/local/intel2013/include
# LIBS =  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lcgns
# F90 = ifort $(ARCH)
# LD  = ifort $(ARCH)

# Gfortran stuff
ARCH = -fautomatic -fmax-stack-var-size=8192 -fconvert=big-endian -fdefault-real-8  -fdefault-double-8 -ffree-line-length-none -L$(LAPACKLIB) -L$(CGNSLIB) -I$(CGNSINCLUDE)
#ARCH = -fautomatic -fmax-stack-var-size=8192  -fdefault-real-8  -fdefault-double-8 -ffree-line-length-none -L$(LAPACKLIB) -L$(CGNSLIB) -I$(CGNSINCLUDE)
LAPACKLIB = $(PWD)/../../lapack
LIBS =  -llapack -lblas -lcgns
LD  = gfortran $(ARCH)
F90 = gfortran $(ARCH)

# CGNS path
CGNSINCLUDE = /usr/local/include
CGNSLIB = /usr/local/lib
# CGNSINCLUDE = $(PWD)/../../libs/include
# CGNSLIB = $(PWD)/../../libs


#=======================================================================
# Standard definitions
OFILES = $(F90FILES:.f90=.o) $(F77FILES:.f=.o) $(CFILES:.c=.o)


VPATH = obj
.SUFFIXES: .out .o .f90 .f .c .inc .h

#=======================================================================
# Targets and dependencies

default:: debug
all:: a.out
debug ::	
	@make all "FLAGS = -fbacktrace  -fcheck=all,no-array-temps $(FLAGS)" # gfortran
	#@make all "FLAGS = -O0 -ftrapuv -gen-interfaces -warn interfaces -heap-arrays -debug all -g -traceback -check bounds -check uninit -check format -check output_conversion -check power -check args $(FLAGS)" # ifort
	
opt::
	@make all "FLAGS =-O2 $(FLAGS)"
static::
	@make all "FLAGS = -static -O2  $(FLAGS)"
clean::
	clear; rm -rf obj a.out *.kmo *.mod *genmod.f90 *genmod.mod $(NAME)

a.out: obj $(OFILES)
	@echo "$(LD) $(OFILES) $(LIBS)"
	@cd obj; $(LD)  $(OFILES) $(LIBS); mv a.out ../$(NAME)
obj:
	mkdir obj

#=======================================================================
# Compilation rules
#
.f90.o:
	$(F90) $(FLAGS) -c $*.f90 -o obj/$*.o
