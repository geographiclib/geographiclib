# Build the test programs...

set (TESTPROGRAMS ProjTest TMTest GeodTest ConicTest NaNTester HarmTest)

# Loop over all the tools, specifying the source and library.
add_custom_target (testprograms ALL)
foreach (TESTPROGRAM ${TESTPROGRAMS})

  add_executable (${TESTPROGRAM} ${TESTPROGRAM}.cpp)
  add_dependencies (testprograms ${TESTPROGRAM})
  target_link_libraries (${TESTPROGRAM} Geographic)

endforeach ()

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
    set_target_properties (${TESTPROGRAMS} PROPERTIES
      INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
  endif ()
endif ()

# Don't install test programs
