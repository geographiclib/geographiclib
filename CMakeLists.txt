# $Id$

cmake_minimum_required (VERSION 2.6)
project (GeographicLib)

set (GeographicLib_VERSION_MAJOR 1)
set (GeographicLib_VERSION_MINOR 14)
set (GeographicLib_VERSION
     "${GeographicLib_VERSION_MAJOR}.${GeographicLib_VERSION_MINOR}")
set (LIBVERSION 9)
set (LIBVERSIONFULL 9.0.1)

if (WIN32)
  option (GEOGRAPHIC_SHARED_LIB "Build GeographicLib as a shared library" OFF)
else (WIN32)
  option (GEOGRAPHIC_SHARED_LIB "Build GeographicLib as a shared library" ON)
endif (WIN32)
option (ENABLE_MATLAB "Try to compile matlab/octave interface" OFF)

if (WIN32)
  set (GEOGRAPHICLIB_GEOID_PATH
    "C:/Documents and Settings/All Users/Application Data/GeographicLib/geoids"
    CACHE STRING "Location for geoid data")
else (WIN32)
  set (GEOGRAPHICLIB_GEOID_PATH
    "${CMAKE_INSTALL_PREFIX}/share/GeographicLib/geoids"
    CACHE STRING "Location for geoid data")
endif (WIN32)

set (CMAKE_DEBUG_POSTFIX _d)

include (CheckTypeSize)
check_type_size ("long double" LONG_DOUBLE BUILTIN_TYPES_ONLY)

include (TestBigEndian)
test_big_endian (WORDS_BIGENDIAN)

configure_file (
    "${PROJECT_SOURCE_DIR}/include/GeographicLib/Config.h.in"
    "${PROJECT_BINARY_DIR}/include/GeographicLib/Config.h" )

if (NOT WIN32)
  find_program (HAVE_POD2MAN pod2man)
  find_program (HAVE_POD2HTML pod2html)
  find_program (HAVE_DOXYGEN doxygen)
endif (NOT WIN32)
if (NOT WIN32 AND HAVE_POD2MAN AND HAVE_POD2HTML AND HAVE_DOXYGEN)
  set (MAINTAINER ON CACHE BOOL "Allow build of additional components")
else (NOT WIN32 AND HAVE_POD2MAN AND HAVE_POD2HTML AND HAVE_DOXYGEN)
  set (MAINTAINER OFF CACHE BOOL "Allow build of additional components")
endif (NOT WIN32 AND HAVE_POD2MAN AND HAVE_POD2HTML AND HAVE_DOXYGEN)

if (ENABLE_MATLAB)
  if (NOT WIN32)
    find_program (HAVE_MEX mkoctfile)
    if (HAVE_MEX)
      set (MEX "mkoctfile")
      set (MEXEXT "mex")
    endif (HAVE_MEX)
  else (NOT WIN32)
    find_program (HAVE_MEX mex.bat)
    if (HAVE_MEX)
      set (MEX "mex")
      set (MEXEXT "mexw32")
    endif (HAVE_MEX)
  endif (NOT WIN32)
else (ENABLE_MATLAB)
  set (HAVE_MEX OFF)
endif (ENABLE_MATLAB)

# Set a default build type for single-configuration
# CMake generators if no build type is set.
if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

if (NOT WIN32)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif (NOT WIN32)

set (TOOLS CartConvert ConicProj GeodesicProj GeoConvert Geod
  GeoidEval MagneticField Planimeter TransverseMercatorProj)

include_directories ("${PROJECT_BINARY_DIR}/include")
include_directories (include)
add_subdirectory (src)
add_subdirectory (include/GeographicLib)
add_subdirectory (tools)
add_subdirectory (man)
add_subdirectory (doc)
add_subdirectory (matlab)
add_subdirectory (python/geographiclib)

if (WIN32)
  set(INSTALL_CMAKE_DIR "cmake")
else (WIN32)
  set(INSTALL_CMAKE_DIR
    "share/cmake/${CMAKE_PROJECT_NAME}-${GeographicLib_VERSION}")
endif (WIN32)

export (TARGETS Geographic ${TOOLS}
  FILE "${PROJECT_BINARY_DIR}/GeographicLibDepends.cmake")

# Don't want cmake to look in build trees, so don't do
# export (PACKAGE GeographicLib)

set (GEOGRAPHICLIB_ROOT_DIR "${PROJECT_BINARY_DIR}")
set (GeographicLib_INCLUDE_DIRS
  "${PROJECT_BINARY_DIR}/include" "${PROJECT_SOURCE_DIR}/include" )
set (GeographicLib_LIBRARY_DIRS "${PROJECT_BINARY_DIR}/src")
set (GEOGRAPHICLIB_CMAKE_DIR "${PROJECT_BINARY_DIR}")
configure_file (cmake/geographiclib-config.cmake.in
  "${PROJECT_BINARY_DIR}/geographiclib-config.cmake" @ONLY)
configure_file (cmake/geographiclib-config-version.cmake.in
  "${PROJECT_BINARY_DIR}/geographiclib-config-version.cmake" @ONLY)

install (EXPORT GeographicLibDepends DESTINATION "${INSTALL_CMAKE_DIR}")

if (WIN32)
  set (GEOGRAPHICLIB_ROOT_DIR "..")
else (WIN32)
  set (GEOGRAPHICLIB_ROOT_DIR "../../..")
endif (WIN32)
set (GEOGRAPHICLIB_CMAKE_DIR "${CMAKE_INSTALL_PREFIX}/${INSTALL_CMAKE_DIR}")
configure_file (cmake/geographiclib-config.cmake.in
  "${PROJECT_BINARY_DIR}/InstallFiles/geographiclib-config.cmake" @ONLY)
configure_file (cmake/geographiclib-config-version.cmake.in
  "${PROJECT_BINARY_DIR}/InstallFiles/geographiclib-config-version.cmake" @ONLY)
install (FILES
  "${PROJECT_BINARY_DIR}/InstallFiles/geographiclib-config.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}")
install (FILES
  "${PROJECT_BINARY_DIR}/InstallFiles/geographiclib-config-version.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}")

if (MAINTAINER)
  add_custom_target (distrib-all)
  add_dependencies (distrib-all distrib-man distrib-doc)
endif (MAINTAINER)

include (InstallRequiredSystemLibraries)
set (CPACK_PACKAGE_VERSION_MAJOR ${GeographicLib_VERSION_MAJOR})
set (CPACK_PACKAGE_VERSION_MINOR ${GeographicLib_VERSION_MINOR})
set (CPACK_PACKAGE_VERSION_PATCH 0)
set (CPACK_PACKAGE_CONTACT charles@karney.com)
set (CPACK_PACKAGE_VENDOR "GeographicLib")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY
  "Geographic library, utilities, and documentation")
set (CPACK_SOURCE_IGNORE_FILES
  "#"
  "~\$"
  "${PROJECT_SOURCE_DIR}/(BUILD|WINDOWS|TAGS|geographiclib-[0-9])"
  "${PROJECT_SOURCE_DIR}/(release-notes|geoids-distrib|geoids-installer)/"
  "${PROJECT_SOURCE_DIR}/(distrib|cgi-bin|scratch|srcL|workspace)/"
  "${PROJECT_SOURCE_DIR}/(tests|testdata)/"
  "${PROJECT_SOURCE_DIR}/(geodesic.*|makefile-admin|install_manifest.txt)"
  "Makefile\$"
  "Config.h\$"
  "TAGS\$"
  "/GeographicLibL"
  "/\\\\.(deps|libs|svn)/"
  "/\\\\.git"
  "/(UpgradeLog|_UpgradeReport)"
  "/(CMakeCache.txt|CPack)"
  "/(CMakeFiles|_CPack_Packages|BACKUP)/"
  "/config\\\\.(h|status)"
  "/Config-ac\\\\.h\$"
  "/(status|libtool)\$"
  "/stamp-"
  "/junk"
  "/octave-core"
  "\\\\.(bak|BAK|patch|tmp|zip|gz|lo|o|ncb|suo|a|la|log|exe|pyc)\$"
  "\\\\.(mex|mexw32|mexa64|kml|kmz)\$"
  "/(Release|Debug|Win32)"
  "\\\\.cache/"
  "proj.*\\\\.user\$"
  "lib.*\\\\.so"
  "tools/[A-Za-z]+\$"
  "tools/geographiclib-get-geoids\$"
  "tools/.*\\\\.usage\$"
  "scripts/geographiclib.js\$"
  "/(fix-permissions|autogen|biblio|js-compress)\\\\.sh"
  "${PROJECT_SOURCE_DIR}/[^/]*\\\\.html" )
set (CPACK_SOURCE_GENERATOR TGZ)
set (CPACK_RESOURCE_FILE_LICENSE ${PROJECT_SOURCE_DIR}/LICENSE.txt)
set (CPACK_PACKAGE_INSTALL_DIRECTORY
  "${CMAKE_PROJECT_NAME}-${GeographicLib_VERSION}")
set (CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY}")

if (WIN32)
  set (CPACK_NSIS_CONTACT "charles@karney.com")
  set (CPACK_NSIS_URL_INFO_ABOUT "http://geographiclib.sf.net")
  set (CPACK_NSIS_HELP_LINK "mailto:charles@karney.com")
  set (CPACK_NSIS_PACKAGE_NAME
    "${CMAKE_PROJECT_NAME} ${GeographicLib_VERSION}")
  set (CPACK_NSIS_DISPLAY_NAME ${CPACK_NSIS_PACKAGE_NAME})
  set (CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY}-win32")
  set (CPACK_NSIS_EXTRA_INSTALL_COMMANDS "
  CreateShortCut \\\"$SMPROGRAMS\\\\$STARTMENU_FOLDER\\\\Library Documentation.lnk\\\" \\\"$INSTDIR\\\\share\\\\doc\\\\GeographicLib\\\\html\\\\index.html\\\"
  CreateShortCut \\\"$SMPROGRAMS\\\\$STARTMENU_FOLDER\\\\Utility Documentation.lnk\\\" \\\"$INSTDIR\\\\share\\\\doc\\\\GeographicLib\\\\html\\\\utilities.html\\\"
")
  set (CPACK_NSIS_EXTRA_UNINSTALL_COMMANDS "
  !insertmacro MUI_STARTMENU_GETFOLDER Application $MUI_TEMP
  Delete \\\"$SMPROGRAMS\\\\$MUI_TEMP\\\\Library Documentation.lnk\\\"
  Delete \\\"$SMPROGRAMS\\\\$MUI_TEMP\\\\Utility Documentation.lnk\\\"
")
  set (CPACK_NSIS_MODIFY_PATH ON)
elseif (APPLE)
  set (CPACK_GENERATOR DMG)
  set (CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY}-darwin")
else (WIN32)
  set (CPACK_GENERATOR TGZ)
endif (WIN32)

include (CPack)

if (MAINTAINER)
  add_custom_target (dist
    COMMAND ${CMAKE_MAKE_PROGRAM} package_source
    COMMAND
      cd _CPack_Packages/Linux-Source/TGZ/${CPACK_SOURCE_PACKAGE_FILE_NAME} &&
      echo include Makefile.mk > Makefile &&
      echo "\\#define" GEOGRAPHICLIB_VERSION_STRING
      \\"${GeographicLib_VERSION}\\" > include/GeographicLib/Config.h
    COMMAND
      cd _CPack_Packages/Linux-Source/TGZ/${CPACK_SOURCE_PACKAGE_FILE_NAME} &&
      touch man/*.usage man/*.1 man/*.1.html && touch doc/html/index.html
    COMMAND
      cd _CPack_Packages/Linux-Source/TGZ &&
      find ${CPACK_SOURCE_PACKAGE_FILE_NAME} -type f |
      tar cfzT ${CMAKE_BINARY_DIR}/${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz -
    COMMAND
      rm -f ${CPACK_SOURCE_PACKAGE_FILE_NAME}.zip &&
      cd _CPack_Packages/Linux-Source/TGZ &&
      find ${CPACK_SOURCE_PACKAGE_FILE_NAME} -type f |
      zip -q ${CMAKE_BINARY_DIR}/${CPACK_SOURCE_PACKAGE_FILE_NAME}.zip -@
  )
  add_dependencies (dist distrib-all)
endif (MAINTAINER)

enable_testing ()
