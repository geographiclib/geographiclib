# Build the test programs...

set (TESTPROGRAMS
  ProjTest TMTest GeodTest ConicTest NaNTester HarmTest EllipticTest GeodShort)

if (NOT WIN32)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
endif ()

# Loop over all the tools, specifying the source and library.
add_custom_target (testprograms ALL)
foreach (TESTPROGRAM ${TESTPROGRAMS})

  add_executable (${TESTPROGRAM} ${TESTPROGRAM}.cpp)
  add_dependencies (testprograms ${TESTPROGRAM})
  target_link_libraries (${TESTPROGRAM} Geographic)

endforeach ()

add_executable (GeodExact GeodExact.cpp
  Geodesic30.cpp GeodesicLine30.cpp
  Geodesic30.hpp GeodesicLine30.hpp)
add_dependencies (testprograms GeodExact)
target_link_libraries (GeodExact Geographic)
set (TESTPROGRAMS ${TESTPROGRAMS} GeodExact)

if (GEOGRAPHIC_SHARED_LIB)
  if (WIN32)
    # Copy the shared library on Windows systems to this directory
    # (tools) so that the tests can be run.
    add_custom_command (TARGET testprograms POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E
	copy $<TARGET_FILE:Geographic> ${CMAKE_CFG_INTDIR}
      COMMENT "Copying shared library to tools directory")
  else ()
    # Set the run time path for shared libraries for non-Windows machines.
    set_target_properties (${TESTPROGRAMS} GeodExact PROPERTIES
      INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
  endif ()
endif ()

# Put all the tools into a folder in the IDE
set_property (TARGET testprograms ${TESTPROGRAMS} PROPERTY FOLDER tests)

# Don't install test programs
