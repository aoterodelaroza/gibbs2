# Suffix for Fortran 90 programs
.SUFFIXES:
.SUFFIXES: .f90 .f $(SUFFIXES)

include ../Makefile.inc

OBJS1=dp1vlu.o dpcoef.o dpolft.o fdump.o i1mach.o j4save.o xercnt.o xerhlt.o\
     xermsg.o xerprn.o xersve.o xgetua.o  

.f90.o:
	$(FC) -c $(FCFLAGS) $<
.f.o:
	$(FC) -c $(FCFLAGS) $<

all: libslatec.a

libslatec.a: $(OBJS1)
	$(AR) rvu libslatec.a $(OBJS1) 

mrproper:
	/bin/rm -f $(OBJS1) libslatec.a 

veryclean:
	/bin/rm -f $(OBJS1) libslatec.a 

clean:
	/bin/rm -f $(OBJS1) 

