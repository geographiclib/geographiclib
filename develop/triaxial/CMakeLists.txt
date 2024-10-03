cmake_minimum_required (VERSION 3.13.0)
project (triaxial)

# -DGeographicLib_DIR=~/geographiclib/BUILD-mpfr -DGEOGRAPHICLIB_PRECISION=5
find_package (GeographicLib 2 REQUIRED)
include_directories (${GeographicLib_INCLUDE_DIRS})
find_package (Boost REQUIRED)

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
set (CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES "$ENV{HOME}/geographiclib/src")

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

add_executable (randgeod randgeod.cpp)
target_link_libraries (randgeod ${GeographicLib_LIBRARIES}
  ${GeographicLib_HIGHPREC_LIBRARIES})

add_executable (TrigfunTest TrigfunTest.cpp)

target_link_libraries (TrigfunTest trilib ${GeographicLib_LIBRARIES}
  ${GeographicLib_HIGHPREC_LIBRARIES})

add_executable (TriaxialTest TriaxialTest.cpp)
target_link_libraries (TriaxialTest trilib ${GeographicLib_LIBRARIES}
  ${GeographicLib_HIGHPREC_LIBRARIES})

add_executable (Geod3Solve Geod3Solve.cpp)
target_link_libraries (Geod3Solve trilib ${GeographicLib_LIBRARIES}
  ${BOOST_LIBRARIES} ${GeographicLib_HIGHPREC_LIBRARIES})

add_executable (TriRand TriRand.cpp)
target_link_libraries (TriRand trilib ${GeographicLib_LIBRARIES}
  ${GeographicLib_HIGHPREC_LIBRARIES})

add_executable (Geod3Test Geod3Test.cpp)
target_link_libraries (Geod3Test trilib ${GeographicLib_LIBRARIES}
  ${BOOST_LIBRARIES} ${GeographicLib_HIGHPREC_LIBRARIES})

add_executable (Geod2Test Geod2Test.cpp)
target_link_libraries (Geod2Test trilib ${GeographicLib_LIBRARIES}
  ${BOOST_LIBRARIES} ${GeographicLib_HIGHPREC_LIBRARIES})
