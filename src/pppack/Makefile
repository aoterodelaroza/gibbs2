# Suffix for Fortran 90 programs
.SUFFIXES:
.SUFFIXES: .f90 .f $(SUFFIXES)

include ../Makefile.inc

AR=ar
OBJS1=cubspl.o interv.o ppvalu.o

.f90.o:
	$(FC) -c $(FCFLAGS) $<
.f.o:
	$(FC) -c $(FCFLAGS) $<

all: libpppack.a

libpppack.a: $(OBJS1)
	$(AR) rvu libpppack.a $(OBJS1) 

mrproper:
	/bin/rm -f $(OBJS1) libpppack.a 

veryclean:
	/bin/rm -f $(OBJS1) libpppack.a 

clean:
	/bin/rm -f $(OBJS1) 


