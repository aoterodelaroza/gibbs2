############# 
# Compiler 
############# 

## The Intel(R) fortran compiler (ifort)
#ifeq ($(DEBUG),1)
#  FC = ifort
#  FCFLAGS = -g -CU -C -traceback -fpe0 -debug
#  LDFLAGS = 
#else
#  FC = ifort
#  FCFLAGS = -O3 
#  LDFLAGS = -O3
#endif

# The GNU fortran compiler (gfortran)
ifeq ($(DEBUG),1)
  FC = gfortran
  FCFLAGS = -g -fbounds-check -Wall -Wunused-parameter -ffpe-trap=invalid -fbacktrace -fdump-core
  LDFLAGS = 
else
  FC = gfortran
  FCFLAGS = -O3
  LDFLAGS = 
endif
