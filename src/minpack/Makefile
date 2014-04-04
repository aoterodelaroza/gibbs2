# Suffix for Fortran 90 programs
.SUFFIXES:
.SUFFIXES: .f90 .f $(SUFFIXES)

include ../Makefile.inc

AR=ar
OBJS1=dpmpar.o enorm.o fdjac2.o lmdif1.o lmdif.o lmpar.o qrfac.o qrsolv.o lmder.o

.f90.o:
	$(FC) -c $(FCFLAGS) $<
.f.o:
	$(FC) -c $(FCFLAGS) $<

all: libminpack.a

libminpack.a: $(OBJS1)
	$(AR) rvu libminpack.a $(OBJS1) 

mrproper:
	/bin/rm -f $(OBJS1) libminpack.a 

veryclean:
	/bin/rm -f $(OBJS1) libminpack.a 

clean:
	/bin/rm -f $(OBJS1) 

