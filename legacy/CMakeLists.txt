cmake_minimum_required (VERSION 2.6)
enable_language (C)
enable_language (Fortran)

# Make the compiler more picky.
if (WIN32)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W4")
else ()
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -pedantic -ansi")
endif ()
add_executable (geodc geodc.c geodesicc.c geodesicc.h)
add_executable (geodf geodf.f geodesicf.f)
if (NOT WIN32)
  target_link_libraries (geodc m)
endif ()

