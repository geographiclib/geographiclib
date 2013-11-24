# Compile (but don't install) a bunch of tiny example programs.  These
# are mainly for including as examples within the doxygen documentation;
# however, compiling them catches most obvious blunders.
file (GLOB EXAMPLE_SOURCES example-*.cpp)
set (EXAMPLES)

foreach (EXAMPLE_SOURCE ${EXAMPLE_SOURCES})
  get_filename_component (EXAMPLE ${EXAMPLE_SOURCE} NAME_WE)
  set (EXAMPLES ${EXAMPLES} ${EXAMPLE})
  add_executable (${EXAMPLE} ${EXAMPLE_SOURCE})
  target_link_libraries (${EXAMPLE} Geographic)
endforeach ()

set (OTHER_EXAMPLES GeoidToGTX.cpp)

foreach (EXAMPLE_SOURCE ${OTHER_EXAMPLES})
  get_filename_component (EXAMPLE ${EXAMPLE_SOURCE} NAME_WE)
  set (EXAMPLES ${EXAMPLES} ${EXAMPLE})
  add_executable (${EXAMPLE} ${EXAMPLE_SOURCE})
  target_link_libraries (${EXAMPLE} Geographic)
endforeach ()

find_package (OpenMP QUIET)

if (OPENMP_FOUND)
  set_target_properties (GeoidToGTX PROPERTIES
    COMPILE_FLAGS ${OpenMP_CXX_FLAGS} COMPILE_DEFINITIONS HAVE_OPENMP=1)
  if (NOT WIN32)
    set_target_properties (GeoidToGTX PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
  endif ()
endif ()

add_custom_target (examples ALL DEPENDS ${EXAMPLES})
add_dependencies (examples Geographic)

if (GEOGRAPHIC_SHARED_LIB AND WIN32)
  # Copy the shared library on Windows systems to this directory (examples)
  # so that the tests can be run.
  add_custom_command (TARGET examples POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E
    copy $<TARGET_FILE:Geographic> ${CMAKE_CFG_INTDIR}
    COMMENT "Copying shared library to examples directory")
endif ()

# Put all the examples into a folder in the IDE
set_property (TARGET examples ${EXAMPLES} PROPERTY FOLDER examples)
