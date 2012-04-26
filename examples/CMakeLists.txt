# $Id$

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

find_package (OpenMP)

if (OPENMP_FOUND OR OpenMP_FOUND)
  set_target_properties (GeoidToGTX PROPERTIES
    COMPILE_FLAGS ${OpenMP_CXX_FLAGS} COMPILE_DEFINITIONS HAVE_OPENMP=1)
  set_target_properties (GeoidToGTX PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
endif (OPENMP_FOUND OR OpenMP_FOUND)

add_custom_target (examples DEPENDS ${EXAMPLES})
add_dependencies (examples Geographic)
