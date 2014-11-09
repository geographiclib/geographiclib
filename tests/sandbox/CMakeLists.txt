cmake_minimum_required (VERSION 2.6)
project (sandbox)

# Set a default build type for single-configuration cmake generators if
# no build type is set.
if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif ()

find_package (GeographicLib 1.25 REQUIRED)

if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories (${GeographicLib_INCLUDE_DIRS})
  add_definitions (${GeographicLib_DEFINITIONS})
endif ()
add_executable (example example-Geodesic.cpp)
target_link_libraries (example ${GeographicLib_LIBRARIES})
