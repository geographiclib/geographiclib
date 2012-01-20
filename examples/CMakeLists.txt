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
endforeach (EXAMPLE_SOURCE)

add_custom_target (examples DEPENDS ${EXAMPLES})
add_dependencies (examples Geographic)
