## enable testing

cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

## list of test bundles
set(TESTDIRS
  001_simple        ## simple inputs and outputs
  002_fits          ## fits using different equations of state
  003_tmodels       ## temperature models
  004_eec           ## empirical energy corrections
  005_elec          ## electronicc contributions to the free energy
  006_phases        ## phase stability
  )

## define the macro that will run the tests
macro (runtests)
  get_filename_component(curdir "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

  foreach(test ${ARGN})
    ## gather info about the test
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cmake")
      include("${CMAKE_CURRENT_SOURCE_DIR}/${test}.cmake")
    else()
      file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/${test}.ing" check_list REGEX "^ *#.*check *: ")
      file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/${test}.ing" delete_list REGEX "^ *#.*delete *: ")
      file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/${test}.ing" label_list REGEX "^ *#.*labels *: ")
    endif()
    string(REGEX REPLACE "^ *#.*labels *:" "" label_list "${label_list}")

    ## copy over the test to the build directory
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${test}.ing DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

    ## gen test
    add_test(
      NAME "GEN:${curdir}/${test}"
      COMMAND gibbs2 -q "${CMAKE_CURRENT_BINARY_DIR}/${test}.ing" "${CMAKE_CURRENT_BINARY_DIR}/${test}.outg"
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
      CONFIGURATIONS gen
      )
    set_property(TEST "GEN:${curdir}/${test}" PROPERTY LABELS "gen ${label_list}")

    ## gen_copy test
    set(_check_list "")
    if (check_list)
      foreach (i IN LISTS check_list)
	string(REGEX REPLACE "^ *#.*check *:" "" i "${i}")
	string(REGEX MATCHALL "([^\ ]+\ |[^\ ]+$)"  llist "${i}")
	list(GET llist 0 file)
	list(APPEND _check_list "${file}")
      endforeach()
      list(TRANSFORM _check_list PREPEND "${CMAKE_CURRENT_BINARY_DIR}/")
      string(REPLACE ";" " " check_string "${_check_list}")

      add_test(
	NAME "GEN-COPY-:${curdir}/${test}"
	COMMAND sh -c "${CMAKE_COMMAND} -E copy ${check_string} ${CMAKE_CURRENT_SOURCE_DIR}/ref/"
	WORKING_DIRECTORY ""
	CONFIGURATIONS gen
	)
      set_property(TEST "GEN-COPY-:${curdir}/${test}" PROPERTY LABELS "gen ${label_list}")
      set_property(TEST "GEN-COPY-:${curdir}/${test}" PROPERTY DEPENDS "GEN:${curdir}/${test}")
    endif()

    ## gen_clean test
    if (_check_list OR delete_list)
      set(_delete_list "")
      foreach (i IN LISTS delete_list)
	string(REGEX REPLACE "^ *#.*delete *:" "" i "${i}")
	string(REGEX MATCHALL "([^\ ]+\ |[^\ ]+$)"  llist "${i}")
	list(APPEND _delete_list "${llist}")
      endforeach()
      list(TRANSFORM _delete_list PREPEND "${CMAKE_CURRENT_BINARY_DIR}/")
      list(APPEND _delete_list "${_check_list}")
      string(REPLACE ";" " " delete_string "${_delete_list}")

      add_test(
	NAME "GEN-CLEAN-:${curdir}/${test}"
	COMMAND sh -c "${CMAKE_COMMAND} -E remove ${delete_string}"
	WORKING_DIRECTORY ""
	CONFIGURATIONS gen
	)
      set_property(TEST "GEN-CLEAN-:${curdir}/${test}" PROPERTY LABELS "gen ${label_list}")
      set_property(TEST "GEN-CLEAN-:${curdir}/${test}" PROPERTY DEPENDS "GEN-COPY:${curdir}/${test}")
    endif()

    ## run test
    add_test(
      NAME "-RUN-:${curdir}/${test}"
      COMMAND gibbs2 -q "${CMAKE_CURRENT_BINARY_DIR}/${test}.ing" "${CMAKE_CURRENT_BINARY_DIR}/${test}.outg"
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
      )
    set_property(TEST "-RUN-:${curdir}/${test}" PROPERTY LABELS "${label_list}")

    ## check tests
    if (NUMDIFF_FOUND AND check_list)
      foreach (i IN LISTS check_list)
	string(REGEX REPLACE "^ *#.*check *:" "" i "${i}")
	string(REGEX MATCHALL "([^\ ]+\ |[^\ ]+$)"  llist "${i}")
	list(LENGTH llist nn)
	list(GET llist 0 file)
	if ("${nn}" GREATER 1)
	  list(GET llist 1 numdiff_opts)
	endif()
	add_test(
	  NAME "CHECK:${curdir}/ref/${file}"
	  COMMAND sh -c "${NUMDIFF_EXE} -s ' \t\n=,:;<>[](){}^' ${numdiff_opts} ${CMAKE_CURRENT_SOURCE_DIR}/ref/${file} ${CMAKE_CURRENT_BINARY_DIR}/${file}"
	  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
	  )
	set_property(TEST "CHECK:${curdir}/ref/${file}" PROPERTY LABELS "${label_list}")
	set_property(TEST "CHECK:${curdir}/ref/${file}" PROPERTY DEPENDS "-RUN-:${curdir}/${test}")
      endforeach()
    endif()
  endforeach()
endmacro()

## create symlink for the test data (build)
if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/dat/")
  message(STATUS "Creating symlink for tests in ${CMAKE_CURRENT_BINARY_DIR}/")
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_CURRENT_SOURCE_DIR}/dat" "${CMAKE_CURRENT_BINARY_DIR}/dat")
endif()

## process all subdirectories
foreach(dir ${TESTDIRS})
  add_subdirectory(${dir})
endforeach()
