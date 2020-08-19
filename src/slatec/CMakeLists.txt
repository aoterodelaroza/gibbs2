set(SOURCES dp1vlu.f dpcoef.f dpolft.f fdump.f i1mach.f j4save.f xercnt.f xerhlt.f xermsg.f
            xerprn.f xersve.f xgetua.f)

add_library(slatec OBJECT ${SOURCES})

