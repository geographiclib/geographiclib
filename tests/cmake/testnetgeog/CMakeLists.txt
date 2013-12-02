cmake_minimum_required (VERSION 2.6)
project (geodesictest)

find_package(GeographicLib 1.34 REQUIRED COMPONENTS NETGeographicLib)

add_executable (${PROJECT_NAME} example-Geodesic-small.cpp)
set_target_properties (${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "/clr")
string (REPLACE "/EHsc" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

get_target_property (NETGeographicLib_LOCATION
  "${GeographicLib_NETGeographicLib_LIBRARIES}" IMPORTED_LOCATION_RELEASE)
set_target_properties (${PROJECT_NAME} PROPERTIES VS_DOTNET_REFERENCES
  "${NETGeographicLib_LOCATION}")

get_target_property (_LIB "${GeographicLib_NETGeographicLib_LIBRARIES}"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE)
get_target_property(_LIBTYPE ${_LIB} TYPE)
if (_LIBTYPE STREQUAL "SHARED_LIBRARY")
  # On Windows systems, copy the shared library to build directory
  if (CMAKE_VERSION VERSION_LESS 2.8.11)
    get_target_property (_LOCR ${_LIB} IMPORTED_LOCATION_RELEASE)
    add_custom_command (TARGET ${PROJECT_NAME} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy "${_LOCR}" ${CMAKE_CFG_INTDIR}
      COMMENT "Copying shared library for GeographicLib")
  else ()
    add_custom_command (TARGET ${PROJECT_NAME} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E
      copy $<TARGET_FILE:${_LIB}> ${CMAKE_CFG_INTDIR}
      COMMENT "Copying shared library for GeographicLib")
  endif ()
endif ()
