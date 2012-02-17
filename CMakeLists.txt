# $Id$

cmake_minimum_required (VERSION 2.6)
project (GeographicLib)

# Version information
set (GeographicLib_VERSION_MAJOR 1)
set (GeographicLib_VERSION_MINOR 18)
set (GeographicLib_VERSION
     "${GeographicLib_VERSION_MAJOR}.${GeographicLib_VERSION_MINOR}")
# The library version tracks the numbering given by libtool in the
# autoconf set up.
set (LIBVERSION 9)
set (LIBVERSIONFULL 9.2.1)

# User-settable cache variables

# (1) Build as a shared library?  On Windows systems, this is typically
# more trouble than it's worth.
if (WIN32)
  option (GEOGRAPHIC_SHARED_LIB "Build GeographicLib as a shared library" OFF)
else (WIN32)
  option (GEOGRAPHIC_SHARED_LIB "Build GeographicLib as a shared library" ON)
endif (WIN32)

# (2) Try to compile the Matlab interfaces?
option (ENABLE_MATLAB "Try to compile matlab/octave interface" OFF)

# (3) Where to look for data files.  Various classes look in the geoids,
# gravity, magnetic, subdirectories of ${GEOGRAPHICLIB_DATA}.
if (WIN32)
  set (GEOGRAPHICLIB_DATA
    "C:/Documents and Settings/All Users/Application Data/GeographicLib"
    CACHE PATH "Location for data for GeographicLib")
else (WIN32)
  set (GEOGRAPHICLIB_DATA
    "${CMAKE_INSTALL_PREFIX}/share/GeographicLib"
    CACHE PATH "Location for data for GeographicLib")
endif (WIN32)

# The debug version of the library is called Geographic_d.
set (CMAKE_DEBUG_POSTFIX _d)

# Determine system properties
include (CheckTypeSize)
check_type_size ("long double" LONG_DOUBLE BUILTIN_TYPES_ONLY)
include (TestBigEndian)
test_big_endian (WORDS_BIGENDIAN)

# Create a Config.h to expose system information to the compiler
configure_file (
    "${PROJECT_SOURCE_DIR}/include/GeographicLib/Config.h.in"
    "${PROJECT_BINARY_DIR}/include/GeographicLib/Config.h" )

# The documentation depends on doxygen.  In addition, the man pages are
# written as pod files and converted to nroff format, C++ code, and
# html.  Because this require tools that may not be available on an
# end-user's system, the creation of the final documentation is
# therefore only done in "MAINTAINER" mode.  The maintainer runs "make
# distrib-all" which installs the transformed documentation files into
# the source tree.
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
if (MAINTAINER)
  add_custom_target (distrib-all)
  add_dependencies (distrib-all distrib-man distrib-doc)
endif (MAINTAINER)

# Look for the tool to compile the Matlab interfaces.  On Linux systems,
# octave is used, since I couldn't find a version of Matlab that
# supported sufficiently recent versions of g++.
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
      execute_process (COMMAND mexext.bat
        OUTPUT_VARIABLE MEXEXT OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif (HAVE_MEX)
  endif (NOT WIN32)
else (ENABLE_MATLAB)
  set (HAVE_MEX OFF)
endif (ENABLE_MATLAB)

# Set a default build type for single-configuration cmake generators if
# no build type is set.
if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

# On non-Windows machine, make the compiler more picky.
if (NOT WIN32)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif (NOT WIN32)

# The list of tools (to be installed into, e.g., /usr/local/bin)
set (TOOLS CartConvert ConicProj GeodesicProj GeoConvert Geod
  GeoidEval Gravity MagneticField Planimeter TransverseMercatorProj)
# The list of scripts (to be installed into, e.g., /usr/local/sbin)
set (SCRIPTS
  geographiclib-get-geoids geographiclib-get-gravity geographiclib-get-magnetic)

# Set the include directories.  Look in ${PROJECT_BINARY_DIR}/include
# first because that's where Config.h will be
include_directories ("${PROJECT_BINARY_DIR}/include")
include_directories (include)

# The list of subdirectories to process
add_subdirectory (src)
add_subdirectory (include/GeographicLib)
add_subdirectory (tools)
add_subdirectory (man)
add_subdirectory (doc)
add_subdirectory (matlab)
add_subdirectory (python/geographiclib)
if (MAINTAINER)
  add_subdirectory (examples)
endif (MAINTAINER)

# config file support for find_package(GeographicLib).  This needs to
# deal with two environments: (1) finding the build tree and (2) finding
# the install tree.  geographiclib-config.cmake detects which situation
# it is handing by looking at @GEOGRAPHICLIB_ROOT_DIR@.  If this is an
# absolute path, it's in the build tree; otherwise, it's in the install
# tree (incidentally, this allows the install tree to be relocated).

# geographiclib-config.cmake for the build tree
set (GEOGRAPHICLIB_ROOT_DIR "${PROJECT_BINARY_DIR}")
set (GeographicLib_INCLUDE_DIRS
  "${PROJECT_BINARY_DIR}/include" "${PROJECT_SOURCE_DIR}/include" )
configure_file (cmake/geographiclib-config.cmake.in
  "${PROJECT_BINARY_DIR}/geographiclib-config.cmake" @ONLY)
configure_file (cmake/geographiclib-config-version.cmake.in
  "${PROJECT_BINARY_DIR}/geographiclib-config-version.cmake" @ONLY)

# geographiclib-config.cmake for the install tree.  It's installed in
# ${INSTALL_CMAKE_DIR} and @GEOGRAPHICLIB_ROOT_DIR@ is the relative
# path to the root from there.
if (WIN32)
  set (INSTALL_CMAKE_DIR "cmake")
  set (GEOGRAPHICLIB_ROOT_DIR "..")
else (WIN32)
  set (INSTALL_CMAKE_DIR
    "share/cmake/${CMAKE_PROJECT_NAME}-${GeographicLib_VERSION}")
  set (GEOGRAPHICLIB_ROOT_DIR "../../..")
endif (WIN32)
# @GeographicLib_INCLUDE_DIRS@ is not used in the install tree;
# reset it to prevent the source and build paths appearing in the
# installed config files
set (GeographicLib_INCLUDE_DIRS IGNORE)
configure_file (cmake/geographiclib-config.cmake.in
  "${PROJECT_BINARY_DIR}/InstallFiles/geographiclib-config.cmake" @ONLY)
configure_file (cmake/geographiclib-config-version.cmake.in
  "${PROJECT_BINARY_DIR}/InstallFiles/geographiclib-config-version.cmake" @ONLY)
install (FILES
  "${PROJECT_BINARY_DIR}/InstallFiles/geographiclib-config.cmake"
  "${PROJECT_BINARY_DIR}/InstallFiles/geographiclib-config-version.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}")

# Make information about the cmake targets (the library and the tools)
# available.
export (TARGETS Geographic ${TOOLS}
  FILE "${PROJECT_BINARY_DIR}/GeographicLibDepends.cmake")
install (EXPORT GeographicLibDepends DESTINATION "${INSTALL_CMAKE_DIR}")

# Don't want find_package to look in build trees, so don't do...
# export (PACKAGE GeographicLib)

# Packaging support; we deal with 
# (1) a source distribution: cmake make a tar.gz file and the zip file
# is created from this.  Only the maintainer can do this, because of
# the need to generate additional documentation files.
# (2) a binary distribution: code is included for Linux, Apple, and
# Windows, but only the Windows distribution has been exercised.

# Need to ensure that system dlls get included in a binary distribution
include (InstallRequiredSystemLibraries)

# The configuration of CPack is via variable that need to be set before
# the include (CPack).
set (CPACK_PACKAGE_VERSION_MAJOR ${GeographicLib_VERSION_MAJOR})
set (CPACK_PACKAGE_VERSION_MINOR ${GeographicLib_VERSION_MINOR})
set (CPACK_PACKAGE_VERSION_PATCH 0)
set (CPACK_PACKAGE_CONTACT charles@karney.com)
set (CPACK_PACKAGE_VENDOR "GeographicLib")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY
  "Geographic library, utilities, and documentation")
# It's necessary to specify a rather long list of files to be excluded
# in the source distribution.
set (CPACK_SOURCE_IGNORE_FILES
  "#"
  "~\$"
  "${PROJECT_SOURCE_DIR}/(BUILD|WINDOWS|TAGS|geographiclib-[0-9])"
  "${PROJECT_SOURCE_DIR}/release-notes/"
  "${PROJECT_SOURCE_DIR}/(geoids-distrib|geoids-installer)/"
  "${PROJECT_SOURCE_DIR}/(magnetic-distrib|magnetic-installer)/"
  "${PROJECT_SOURCE_DIR}/(gravity-distrib|gravity-installer)/"
  "${PROJECT_SOURCE_DIR}/(magnetic|gravity)/"
  "${PROJECT_SOURCE_DIR}/(distrib|cgi-bin|scratch|srcL|workspace)/"
  "${PROJECT_SOURCE_DIR}/(tests|testdata)/"
  "${PROJECT_SOURCE_DIR}/python/(MANIFEST\$|build/)"
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
  "\\\\.(mex|mexw32|mexw64|mexa64|kml|kmz)\$"
  "/(Release|Debug|Win32)"
  "\\\\.cache/"
  "proj.*\\\\.user\$"
  "lib.*\\\\.so"
  "tools/[A-Za-z]+\$"
  "tools/geographiclib-get-geoids\$"
  "tools/geographiclib-get-magnetic\$"
  "tools/geographiclib-get-gravity\$"
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
  # The Windows binary packager is NSIS.  Set the necessary variables
  # for this.
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
  # Not tested
  set (CPACK_GENERATOR DMG)
  set (CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY}-darwin")
else (WIN32)
  # Not tested
  set (CPACK_GENERATOR TGZ)
endif (WIN32)

include (CPack)

# Another maintainer-specific target is building the source distribution
# via the target dist.  This calls package_source to make a tar.gz file.
# However this needs to be touched up to support the vanilla Makefiles
# provided with GeographicLib.  This entails
# (1) creating Makefile (which includes Makefile.mk);
# (2) creating a bare-bones Config.h (with just the version information);
# (3) making sure that make thinks the generated documentation files are
# up-to-date.
# Then a new tar.gz file and zip file are created.  To avoid potential
# problems with directory permissions, tar and zip are told only to
# archive the files.
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

# Add a test target; the tests are in tools.
enable_testing ()
