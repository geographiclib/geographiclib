cmake_minimum_required (VERSION 3.13.0)
project (triaxial)

# Adding files from repo triaxial/src to geographiclib/src/triaxial
# This follows https://stackoverflow.com/a/11426261/837055

# cd ~/git/triaxial
# git log --pretty=email --patch-with-stat --reverse --full-index --binary -m --first-parent -- src > /tmp/patch

# cd /tmp
# git clone ~/geographiclib
# cd geographiclib
# git am --committer-date-is-author-date -p2 --directory=src/triaxial < /tmp/patch

# Check everything is OK.
# git push

# -DGeographicLib_DIR=~/geographiclib/BUILD-mpfr -DGEOGRAPHICLIB_PRECISION=5
find_package (GeographicLib 2.6 CONFIG REQUIRED)
find_package (Boost 1.70 CONFIG REQUIRED)

if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  # Set a default build type for single-configuration cmake generators
  # if no build type is set.
  set (CMAKE_BUILD_TYPE "Release")
endif ()

if (NOT MSVC)
  set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif ()

set (CONVERT_WARNINGS_TO_ERRORS ON)
# We require C++17
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

# For kissfft.hh
if (APPLE)
  # clang needs Boost include directories too
  set (CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES
    "${GeographicLib_INCLUDE_DIRS}/../src" "${Boost_INCLUDE_DIRS}")
else ()
  set (CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES
    "${GeographicLib_INCLUDE_DIRS}/../src")
endif ()

# Make the compiler more picky.
include (CheckCXXCompilerFlag)
if (MSVC)
  string (REGEX REPLACE "/W[0-4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  # Turn on parallel builds for Visual Studio
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /std:c++latest /MP")
else ()
  # Don't treat creation of NaNs as exceptions.
  # N.B., flags -fno-rounding-math -fno-signaling-nans included by default
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-trapping-math -fno-math-errno")
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

set (BUILD_SHARED_LIBS ON)

set (SOURCES
  Trigfun.cpp
  Triaxial.cpp
  TriaxialLine.cpp
  Angle.cpp
  TriaxialODE.cpp)

set (HEADERS
  Trigfun.hpp
  Triaxial.hpp
  TriaxialLine.hpp
  Angle.hpp
  TriaxialODE.hpp)

add_library (trilib ${SOURCES} ${HEADERS})
target_link_libraries (trilib ${GeographicLib_LIBRARIES}
  ${GeographicLib_HIGHPREC_LIBRARIES} ${BOOST_LIBRARIES})

add_executable (randgeod randgeod.cpp)
target_link_libraries (randgeod trilib)

add_executable (TrigfunTest TrigfunTest.cpp)
target_link_libraries (TrigfunTest trilib)

add_executable (TriaxialTest TriaxialTest.cpp)
target_link_libraries (TriaxialTest trilib)

add_executable (Geod3Solve Geod3Solve.cpp)
target_link_libraries (Geod3Solve trilib)

add_executable (Geod3ODE Geod3ODE.cpp)
target_link_libraries (Geod3ODE trilib)

add_executable (TriRand TriRand.cpp)
target_link_libraries (TriRand trilib)

add_executable (Geod3Test Geod3Test.cpp)
target_link_libraries (Geod3Test trilib)

add_executable (InvDump InvDump.cpp)
target_link_libraries (InvDump trilib)

if (GEOGRAPHICLIB_PRECISION LESS 5)
  add_executable (odetest odetest.cpp)
  target_link_libraries (odetest
    ${BOOST_LIBRARIES} ${GeographicLib_HIGHPREC_LIBRARIES})
endif ()
