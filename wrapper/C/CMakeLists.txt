project (geoidtest)

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
  string (REGEX REPLACE "/W[0-4]" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W4")
else ()
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
  # check for C++11 support.
  include (CheckCXXCompilerFlag)
  foreach (_F 1z 14 1y 11 0x)
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

find_package (GeographicLib REQUIRED COMPONENTS SHARED)
include_directories (${GeographicLib_INCLUDE_DIRS})

add_executable (${PROJECT_NAME} ${PROJECT_NAME}.c cgeoid.cpp)
target_link_libraries (${PROJECT_NAME} ${GeographicLib_LIBRARIES})

get_target_property (GEOGRAPHICLIB_LIB_TYPE ${GeographicLib_LIBRARIES} TYPE)
if (GEOGRAPHICLIB_LIB_TYPE STREQUAL "SHARED_LIBRARY")
  if (WIN32)
    add_custom_command (TARGET ${PROJECT_NAME} POST_BUILD
      COMMAND
        ${CMAKE_COMMAND} -E
        copy $<TARGET_FILE:${GeographicLib_LIBRARIES}> ${CMAKE_CFG_INTDIR}
      COMMENT "Installing shared libraries in build tree")
  else ()
    # Set the run time path for shared libraries for non-Windows machines.
    set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif ()
endif ()
