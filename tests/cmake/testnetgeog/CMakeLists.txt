project (geodesictest)
cmake_minimum_required (VERSION 2.8.7) # required for VS_DOTNET_REFERENCES

find_package (GeographicLib 1.34 REQUIRED COMPONENTS NETGeographicLib)

add_executable (${PROJECT_NAME} example-Geodesic-small.cpp)
set_target_properties (${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "/clr")
string (REPLACE "/RTC1" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
string (REPLACE "/EHsc" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

# This is set up for Release builds only.  Change RELEASE to DEBUG for
# Debug builds.
get_target_property (_LOC "${GeographicLib_NETGeographicLib_LIBRARIES}"
  IMPORTED_LOCATION_RELEASE)
set_target_properties (${PROJECT_NAME} PROPERTIES VS_DOTNET_REFERENCES ${_LOC})

get_target_property (_LIB "${GeographicLib_NETGeographicLib_LIBRARIES}"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE)
get_target_property (_LIBTYPE ${_LIB} TYPE)
if (_LIBTYPE STREQUAL "SHARED_LIBRARY")
  # On Windows systems, copy the shared library to build directory
  add_custom_command (TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E
    copy $<TARGET_FILE:${_LIB}> ${CMAKE_CFG_INTDIR}
    COMMENT "Copying shared library for GeographicLib")
endif ()
