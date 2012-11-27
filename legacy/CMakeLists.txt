cmake_minimum_required (VERSION 2.6)

if (WIN32)
  project (GeographicLib-legacy C)
else ()
  project (GeographicLib-legacy C Fortran)
endif ()

# Set a default build type for single-configuration cmake generators if
# no build type is set.
if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif ()

# Make the compiler more picky.
if (WIN32)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W4")
else ()
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -pedantic -ansi")
endif ()
add_executable (geodc geodc.c geodesicc.c geodesicc.h)
if (NOT WIN32)
  target_link_libraries (geodc m)
endif ()

if (NOT WIN32)
  set (CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS} -Wall -Wextra -pedantic -std=f95 -fimplicit-none")
  add_executable (geodf geodf.f geodesicf.f)
endif ()

