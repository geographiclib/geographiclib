cmake_minimum_required (VERSION 3.1.0)
project (sandbox)

# Set a default build type for single-configuration cmake generators if
# no build type is set.
if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif ()

find_package (GeographicLib 1.39 REQUIRED)

add_executable (example example-Geodesic.cpp)
target_link_libraries (example ${GeographicLib_LIBRARIES})
