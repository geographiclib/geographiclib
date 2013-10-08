# Compile (but don't install) a bunch of tiny example programs for
# NETGeographic.  These are mainly for including as examples within the
# doxygen documentation; however, compiling them catches most obvious
# blunders.

file (GLOB EXAMPLE_SOURCES example-*.cpp)
set (EXAMPLES)

foreach (EXAMPLE_SOURCE ${EXAMPLE_SOURCES})
  get_filename_component (EXAMPLE ${EXAMPLE_SOURCE} NAME_WE)
  set (EXAMPLE "net${EXAMPLE}")
  set (EXAMPLES ${EXAMPLES} ${EXAMPLE})
  add_executable (${EXAMPLE} ${EXAMPLE_SOURCE})
  add_dependencies (${EXAMPLE} NETGeographic)
  set_target_properties (${EXAMPLE} PROPERTIES COMPILE_FLAGS "/clr")
  set_target_properties (${EXAMPLE} PROPERTIES VS_DOTNET_REFERENCES
    "${CMAKE_CURRENT_BINARY_DIR}/../../NETGeographicLib/${CMAKE_CFG_INTDIR}/NETGeographic.dll")
endforeach ()

if (CMAKE_CXX_FLAGS_DEBUG MATCHES "/RTC1")
  string (REPLACE "/RTC1" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
endif ()

if (CMAKE_CXX_FLAGS MATCHES "/EHsc")
  string(REPLACE "/EHsc" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif ()

add_custom_target (netexamples ALL DEPENDS ${EXAMPLES})
add_dependencies (netexamples NETGeographic)

if (GEOGRAPHIC_SHARED_LIB AND WIN32)
  # Copy the shared library on Windows systems to this directory (netexamples)
  # so that the tests can be run.
  add_custom_command (TARGET netexamples POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E
    copy $<TARGET_FILE:Geographic> ${CMAKE_CFG_INTDIR}
    COMMENT "Copying shared library to netexamples directory")
endif ()

# Put all the examples into a folder in the IDE
set_property (TARGET netexamples ${EXAMPLES} PROPERTY FOLDER .NET-examples)
