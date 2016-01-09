project (PyGeographicLib)

cmake_minimum_required (VERSION 2.8.8)

# Set a default build type for single-configuration cmake generators if
# no build type is set.
if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif ()

# Make the compiler more picky.
if (MSVC)
  string (REGEX REPLACE "/W[0-4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
else ()
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
  # check for C++11 support.
  include (CheckCXXCompilerFlag)
  foreach (_F 11 0x)
    set (CXX11_FLAG "-std=c++${_F}")
    set (_T CXX11TEST${_F})
    check_cxx_compiler_flag (${CXX11_FLAG} ${_T})
    if (${_T})
      break ()
    endif ()
    unset (CXX11_FLAG)
  endforeach ()
  if (CXX11_FLAG)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX11_FLAG}")
  endif ()
endif ()

# The version of python that boost-python uses.  Also used for the
# installation directory.
set (PYTHON_VERSION 2.7)

# Requires python + python devel
find_package (PythonLibs ${PYTHON_VERSION} REQUIRED)

# Required boost-python + boost-devel
find_package (Boost REQUIRED COMPONENTS python)

find_package (GeographicLib REQUIRED COMPONENTS SHARED)

include_directories (${GeographicLib_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_library (${PROJECT_NAME} MODULE ${PROJECT_NAME}.cpp)

# Don't include the "lib" prefix on the output name
set_target_properties (${PROJECT_NAME} PROPERTIES PREFIX "")
target_link_libraries (${PROJECT_NAME} ${GeographicLib_LIBRARIES}
  ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

install (TARGETS ${PROJECT_NAME} LIBRARY
  # if CMAKE_INSTALL_PREFIX=~/.local then this specifies a directory in
  # the default path.
  DESTINATION lib/python${PYTHON_VERSION}/site-packages)