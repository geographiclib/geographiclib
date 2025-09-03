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
find_package (Boost 1.70 CONFIG)
if (Boost_FOUND)
  add_definitions (-DHAVE_BOOST=1)
  include_directories ("${Boost_INCLUDE_DIRS}")
endif ()

if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  # Set a default build type for single-configuration cmake generators
  # if no build type is set.
  set (CMAKE_BUILD_TYPE "Release")
endif ()

if (NOT MSVC)
  set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif ()

set (CONVERT_WARNINGS_TO_ERRORS ON)
set (REQUIRED_COMPILER_FEATURES "cxx_std_17")

# For kissfft.hh
set (CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES
  "${GeographicLib_INCLUDE_DIRS}/../src")

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

set (SOURCES
  Trigfun.cpp
  Triaxial.cpp
  TriaxialLine.cpp
  TriaxialCartesian.cpp
  Angle.cpp
  TriaxialODE.cpp)

set (HEADERS
  Trigfun.hpp
  Triaxial.hpp
  TriaxialLine.hpp
  TriaxialCartesian.hpp
  Angle.hpp
  TriaxialODE.hpp)

if (MSVC)
  set (BUILD_SHARED_LIBS OFF)
else ()
  set (BUILD_SHARED_LIBS ON)
endif ()

add_library (trilib ${SOURCES} ${HEADERS})
target_link_libraries (trilib ${GeographicLib_LIBRARIES}
  ${GeographicLib_HIGHPREC_LIBRARIES} ${BOOST_LIBRARIES})
target_compile_features (trilib PUBLIC ${REQUIRED_COMPILER_FEATURES})

add_executable (TrigfunTest TrigfunTest.cpp)
target_link_libraries (TrigfunTest trilib)

add_executable (Geod3Solve Geod3Solve.cpp)
target_link_libraries (Geod3Solve trilib)

add_executable (Cart3Convert Cart3Convert.cpp)
target_link_libraries (Cart3Convert trilib)

add_executable (Geod3Test Geod3Test.cpp)
target_link_libraries (Geod3Test trilib)

if (Boost_FOUND)
  add_executable (Geod3ODE Geod3ODE.cpp)
  target_link_libraries (Geod3ODE trilib)
endif ()
