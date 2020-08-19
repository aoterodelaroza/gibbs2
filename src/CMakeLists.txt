## sources
set(SOURCES_Fortran debye.f90 dgam.f eos.f90 evfunc.f90 fit.f90 gnuplot_templates.f90 
                    param.f90 setvariables.f90 tools.f90 topcalc.f90 varbas.f90)
set(SOURCES_Fortran_Main gibbs2.f90)
set(SOURCES ${SOURCES_Fortran} ${SOURCES_Fortran_Main})

## executable
add_executable(gibbs2 ${SOURCES})

## minpack
add_subdirectory(minpack)
include_directories("${CMAKE_BINARY_DIR}/src/minpack")
target_sources(gibbs2 PRIVATE $<TARGET_OBJECTS:minpack>)

## pppack
add_subdirectory(pppack)
include_directories("${CMAKE_BINARY_DIR}/src/pppack")
target_sources(gibbs2 PRIVATE $<TARGET_OBJECTS:pppack>)

## slatec
add_subdirectory(slatec)
include_directories("${CMAKE_BINARY_DIR}/src/slatec")
target_sources(gibbs2 PRIVATE $<TARGET_OBJECTS:slatec>)

## install the binary
install(TARGETS gibbs2 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
