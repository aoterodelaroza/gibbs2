set(SOURCES dpmpar.f enorm.f fdjac2.f lmder.f lmdif1.f lmdif.f lmpar.f qrfac.f qrsolv.f)
add_library(minpack OBJECT ${SOURCES})

