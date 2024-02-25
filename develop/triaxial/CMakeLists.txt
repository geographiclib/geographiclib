project (trixial)
cmake_minimum_required (VERSION 3.13.0)

# -DGeographicLib_DIR=~/geographiclib/BUILD-mpfr -DGEOGRAPHICLIB_PRECISION=5
find_package (GeographicLib 2 REQUIRED)

if (NOT MSVC)
  set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif ()

set (CONVERT_WARNINGS_TO_ERRORS ON)

# Make the compiler more picky.
include (CheckCXXCompilerFlag)
if (MSVC)
  string (REGEX REPLACE "/W[0-4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  # Turn on parallel builds for Visual Studio
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /std:c++latest /MP")
else ()
  set (FLOAT_CONVERSION_FLAG "-Wfloat-conversion")
  check_cxx_compiler_flag (${FLOAT_CONVERSION_FLAG} FLOAT_CONVERSION)
  if (NOT FLOAT_CONVERSION)
    set (FLOAT_CONVERSION_FLAG)
  endif ()
  set (CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -Wall -Wextra ${FLOAT_CONVERSION_FLAG} -Wunused")
endif ()
if (CONVERT_WARNINGS_TO_ERRORS)
  if (MSVC)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /WX")
    set (CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} /WX")
    set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /WX")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /WX")
  else ()
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
  endif ()
endif ()

add_executable (randgeod randgeod.cpp)
target_link_libraries (randgeod ${GeographicLib_LIBRARIES}
  ${GeographicLib_HIGHPREC_LIBRARIES})



