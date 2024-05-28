## project and version
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

## some project variables
file(READ "VERSION" gibbs2_VERSION)
string(REGEX REPLACE " *\n$" "" gibbs2_VERSION "${gibbs2_VERSION}")
set(gibbs2_URL "https://aoterodelaroza.github.io/gibbs2/")
set(gibbs2_EMAIL "aoterodelaroza@gmail.com")
set(gibbs2_LICENSE "GNU/GPL version 3")
set(gibbs2_DESCRIPTION "Thermodynamics of solids under pressure and temperature in the quasiharmonic approximation")
string(TIMESTAMP DATE "%a %d %b %Y %H:%M:%S")

## project definition
project(gibbs2 VERSION "${gibbs2_VERSION}" LANGUAGES Fortran)

## find git to see if we can have a git commit; n/a otherwise
find_package(Git QUIET)
if (GIT_FOUND)
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" describe --always HEAD
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    OUTPUT_VARIABLE gibbs2_GIT_COMMIT
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()
if (NOT gibbs2_GIT_COMMIT)
  set(gibbs2_GIT_COMMIT "n/a")
endif()

## global flags for cmake
set(CMAKE_VERBOSE_MAKEFILE ON)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)
include(GNUInstallDirs)
include(utils)

## default build type, and type options
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release)

## interface options
option(BUILD_STATIC "Build gibbs2 statically." OFF)

## testing
option(BUILD_TESTING "Enable the regression tests for the gibbs2 build." OFF)
include(CTest)

## perform a static build
if (BUILD_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
  find_library(LIBDL_LIBRARY NAMES dl)
  if(NOT LIBDL_LIBRARY)
    error("Dynamic linking static library not found")
  endif()
  mark_as_advanced(LIBDL_LIBRARY)
  set(CMAKE_Fortran_LINK_EXECUTABLE "${CMAKE_Fortran_LINK_EXECUTABLE} ${LIBDL_LIBRARY}")
endif()

## compiler flags, fortran, and check versions
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  set(FFLAGS_DEBUG "-g -O0 -C -traceback -debug all -fpp -check all")
  set(FFLAGS_RELEASE "-O3")
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  set(FFLAGS_DEBUG "-g -Og -fcheck=all -fbounds-check -Wall -Wunused-parameter -Wno-maybe-uninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace -fdump-core -cpp")
  set(FFLAGS_RELEASE "-O3")
endif()

## set the compiler flags
string(TOUPPER ${CMAKE_BUILD_TYPE} BT)
set(CMAKE_Fortran_FLAGS_${BT} "${FFLAGS_${BT}}")

## numdiff
find_package(NUMDIFF)

## the syntax.txt target
add_custom_target(syntax 
  COMMAND wget -qO - 'raw.githubusercontent.com/aoterodelaroza/aoterodelaroza.github.io/master/_gibbs2/12_syntax.md' | awk '/~~~/{a=!a\;next}a{print}' > syntax.txt
  )

## process subdirectories
add_subdirectory(src)
if (BUILD_TESTING AND NUMDIFF_FOUND)
  add_subdirectory(tests)
endif()

## summary
set(_message
"\n\n### ${CMAKE_PROJECT_NAME}, ${CMAKE_PROJECT_VERSION} (git:${gibbs2_GIT_COMMIT})###\n"
"System: ${CMAKE_SYSTEM}\n"
"Date: ${DATE}\n"
"Build type: ${CMAKE_BUILD_TYPE}\n"
"Fortran compiler: ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID}, ${CMAKE_Fortran_COMPILER_VERSION})\n"
"Fortran flags: ${CMAKE_Fortran_FLAGS_${BT}}\n"
"Install directory: ${CMAKE_INSTALL_PREFIX}\n"
)
if (NUMDIFF_FOUND)
  set(_message "${_message}Numdiff: ${NUMDIFF_EXE}\n")
else()
  set(_message "${_message}Numdiff: not found\n")
endif()
set(_message "${_message}################################\n")
message(STATUS ${_message})
