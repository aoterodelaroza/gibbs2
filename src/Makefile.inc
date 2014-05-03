############# 
# Compiler 
############# 

# The Intel(R) fortran compiler (ifort)
ifeq ($(DEBUG),1)
  FC = ifort
  FCFLAGS = -g -CU -C -traceback -fpe0 -debug
  LDFLAGS = 
  AR = ar
  EXE =
else
  FC = ifort
  FCFLAGS = -O3 
  LDFLAGS = -O3
  AR = ar
  EXE =
endif

## The GNU fortran compiler (gfortran)
#ifeq ($(DEBUG),1)
#  FC = gfortran
#  FCFLAGS = -g -fbounds-check -Wall -Wunused-parameter -ffpe-trap=invalid -fbacktrace -fdump-core
#  LDFLAGS = 
#  AR = ar
#  EXE =
#else
#  FC = gfortran
#  FCFLAGS = -O3
#  LDFLAGS = 
#  AR = ar
#  EXE =
#endif

## The minimalistic GNU cross-compiler (windows target)
#ifeq ($(DEBUG),1)
#  FC = i686-w64-mingw32-gfortran
#  FCFLAGS = -g -fbounds-check -Wall -Wunused-parameter -ffpe-trap=invalid -fbacktrace -fdump-core
#  LDFLAGS = -static-libgfortran -static-libgcc
#  AR = i686-w64-mingw32-ar
#  EXE = .exe
#else
#  FC = i686-w64-mingw32-gfortran
#  FCFLAGS = -O3
#  LDFLAGS = -static-libgfortran -static-libgcc
#  AR = i686-w64-mingw32-ar
#  EXE = .exe
#endif
