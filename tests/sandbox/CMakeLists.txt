cmake_minimum_required (VERSION 3.7.0)
project (sandbox)

# Set a default build type for single-configuration cmake generators if
# no build type is set.
if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif ()

find_package (GeographicLib 2.0 REQUIRED)
message(STATUS "GeographicLib_LIBRARIES: ${GeographicLib_LIBRARIES}")

add_executable (sandbox sandbox.cpp)
target_link_libraries (sandbox ${GeographicLib_LIBRARIES})
