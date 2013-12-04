# Build the test programs...

set (TESTPROGRAMS
  ProjTest TMTest GeodTest ConicTest NaNTester HarmTest EllipticTest)

if (NOT WIN32 AND CXX11FLAG)
  set (TESTPROGRAMS ${TESTPROGRAMS} GeodShort)
endif ()

# Loop over all the tools, specifying the source and library.
add_custom_target (testprograms)
foreach (TESTPROGRAM ${TESTPROGRAMS})

  add_executable (${TESTPROGRAM} EXCLUDE_FROM_ALL ${TESTPROGRAM}.cpp)
  add_dependencies (testprograms ${TESTPROGRAM})
  target_link_libraries (${TESTPROGRAM} ${PROJECT_LIBRARIES})

endforeach ()

add_executable (GeodExact EXCLUDE_FROM_ALL GeodExact.cpp
  Geodesic30.cpp GeodesicLine30.cpp
  Geodesic30.hpp GeodesicLine30.hpp)
add_dependencies (testprograms GeodExact)
target_link_libraries (GeodExact ${PROJECT_LIBRARIES})
set (TESTPROGRAMS ${TESTPROGRAMS} GeodExact)

if (MSVC)
  get_target_property (_LIBTYPE ${PROJECT_LIBRARIES} TYPE)
  if (_LIBTYPE STREQUAL "SHARED_LIBRARY")
    # Copy the shared library on Windows systems to this directory
    # (examples) so that the tests can be run.
    if (CMAKE_VERSION VERSION_LESS 2.8.11)
      get_target_property (_LOCR ${PROJECT_LIBRARIES} LOCATION_RELEASE)
      get_target_property (_LOCD ${PROJECT_LIBRARIES} LOCATION_DEBUG)
      add_custom_command (TARGET examples POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy "${_LOCR}" ${CMAKE_CFG_INTDIR}
        COMMAND ${CMAKE_COMMAND} -E copy "${_LOCD}" ${CMAKE_CFG_INTDIR}
        COMMENT "Copying shared library to examples directory")
    else ()
      add_custom_command (TARGET examples POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E
        copy $<TARGET_FILE:${PROJECT_LIBRARIES}> ${CMAKE_CFG_INTDIR}
        COMMENT "Copying shared library to examples directory")
    endif ()
  endif ()
endif ()

# Put all the tools into a folder in the IDE
set_property (TARGET testprograms ${TESTPROGRAMS} PROPERTY FOLDER tests)

# Don't install test programs
