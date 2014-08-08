cmake_minimum_required (VERSION 2.6)

project (GeographicLib-legacy-Fortran Fortran)

# Set a default build type for single-configuration cmake generators if
# no build type is set.
if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif ()

set (CMAKE_Fortran_FLAGS
  "${CMAKE_Fortran_FLAGS} -Wall -Wextra -pedantic -std=f95 -fimplicit-none")
# -Wextra now turns on warnings about comparing reals for equality
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-compare-reals")

set (TOOLS geoddirect geodinverse planimeter)

foreach (TOOL ${TOOLS})
  add_executable (${TOOL} ${TOOL}.for geodesic.for geodesic.inc)
endforeach ()
