cmake_minimum_required (VERSION 3.13.0)
project (trixial)

# -DGeographicLib_DIR=~/geographiclib/BUILD-mpfr -DGEOGRAPHICLIB_PRECISION=5
find_package (GeographicLib 2 REQUIRED)
include_directories (${GeographicLib_INCLUDE_DIRS})
find_package (Boost REQUIRED)

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

add_executable (randgeod randgeod.cpp)
target_link_libraries (randgeod ${GeographicLib_LIBRARIES}
  ${GeographicLib_HIGHPREC_LIBRARIES})

add_executable (TrigfunTest TrigfunTest.cpp Trigfun.cpp Trigfun.hpp
   Triaxial.cpp Triaxial.hpp)
target_link_libraries (TrigfunTest ${BOOST_LIBRARIES} ${GeographicLib_LIBRARIES}
  ${GeographicLib_HIGHPREC_LIBRARIES})

add_executable (TriaxialTest TriaxialTest.cpp Trigfun.cpp Trigfun.hpp
   Triaxial.cpp Triaxial.hpp TriaxialLine.cpp TriaxialLine.hpp)
target_link_libraries (TriaxialTest ${BOOST_LIBRARIES} ${GeographicLib_LIBRARIES}
  ${GeographicLib_HIGHPREC_LIBRARIES})

add_executable (Geod3Solve Geod3Solve.cpp Trigfun.cpp Trigfun.hpp
   Triaxial.cpp Triaxial.hpp TriaxialLine.cpp TriaxialLine.hpp)
target_link_libraries (Geod3Solve ${BOOST_LIBRARIES} ${GeographicLib_LIBRARIES}
  ${GeographicLib_HIGHPREC_LIBRARIES})
