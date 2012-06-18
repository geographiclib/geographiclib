cmake_minimum_required (VERSION 2.6)
project (GeographicLib)

# Version information
set (GeographicLib_VERSION_MAJOR 1)
set (GeographicLib_VERSION_MINOR 23)
set (GeographicLib_VERSION
  "${GeographicLib_VERSION_MAJOR}.${GeographicLib_VERSION_MINOR}")
set (PACKAGE_VERSION ${GeographicLib_VERSION})
# The library version tracks the numbering given by libtool in the
# autoconf set up.
set (LIBVERSION 9)
set (LIBVERSIONFULL 9.5.0)

# User-settable variables

# (1) COMMON_INSTALL_PATH governs the installation convention.  If it
# is on ON (the Linux default), the installation is to a common
# directory, e.g., /usr/local.  If it is OFF (the Windows default),
# the installation directory contains the package name, e.g.,
# C:/pkg/GeographicLib-1.22.  The installation directories for the
# documentation, cmake configuration, python and matlab interfaces all
# depend on the variable with deeper paths relative to
# CMAKE_INSTALL_PREFIX being used when it's ON.

if (WIN32)
  option (COMMON_INSTALL_PATH "Use a common installation path for packages" OFF)
else ()
  option (COMMON_INSTALL_PATH "Use a common installation path for packages" ON)
endif ()

# (2) PACKAGE_PATH and INSTALL_PATH govern the find_package search
# path and the installation directory.  (find_package is not used by
# GeographicLib since it doesn't depend on other packages.  However
# PACKAGE_PATH is used here for uniformity with other packages which
# adopt the same conventions.)
#
# If PACKAGE_PATH is defined, it is prepended to CMAKE_PREFIX_PATH.
#
# If INSTALL_PATH is not specified but PACKAGE_PATH is, then
# INSTALL_PATH is set to
#   ${PACKAGE_PATH}, if COMMON_INSTALL_PATH is ON;
#   ${PACKAGE_PATH}/${PROJECT_NAME}-${PACKAGE_VERSION}, otherwise.
#
# If INSTALL_PATH is now defined, then set CMAKE_INSTALL_PREFIX to
# INSTALL_PATH.
#
# Typically, only PACKAGE_PATH needs to be specified, e.g.,
# cmake -D PACKAGE_PATH=/opt .. (on Linux)
#   => CMAKE_PREFIX_PATH=/opt   CMAKE_INSTALL_PREFIX=/opt
# cmake -D PACKAGE_PATH=C:/pkg .. (on Windows)
#   => CMAKE_PREFIX_PATH=C:/pkg CMAKE_INSTALL_PREFIX=C:/pkg/GeographicLib-1.22

if (DEFINED PACKAGE_PATH)
  set (CMAKE_PREFIX_PATH ${PACKAGE_PATH} ${CMAKE_PREFIX_PATH})
  message (STATUS "CMAKE_PREFIX_PATH set to ${CMAKE_PREFIX_PATH}")
elseif (DEFINED CMAKE_PREFIX_PATH
    AND CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  # Backwards compatibility (support the mechanism advertised in
  # versions 1.20 and 1.21).
  list (GET CMAKE_PREFIX_PATH 0 PACKAGE_PATH)
endif ()

if (NOT DEFINED INSTALL_PATH AND DEFINED PACKAGE_PATH)
  if (COMMON_INSTALL_PATH)
    set (INSTALL_PATH ${PACKAGE_PATH} CACHE PATH "Installation directory" FORCE)
  else ()
    set (INSTALL_PATH ${PACKAGE_PATH}/${PROJECT_NAME}-${PACKAGE_VERSION}
      CACHE PATH "Installation directory" FORCE)
  endif ()
endif ()
if (DEFINED INSTALL_PATH)
  set (CMAKE_INSTALL_PREFIX ${INSTALL_PATH})
endif ()
message (STATUS "CMAKE_INSTALL_PREFIX set to ${CMAKE_INSTALL_PREFIX}")

# (3) Build as a shared library?  On Windows systems, this is typically
# more trouble than it's worth.
if (WIN32)
  option (GEOGRAPHIC_SHARED_LIB "Build GeographicLib as a shared library" OFF)
else ()
  option (GEOGRAPHIC_SHARED_LIB "Build GeographicLib as a shared library" ON)
endif ()

# (4) Compile the Matlab interfaces?  Skip Matlab compilation if OFF
set (MATLAB_COMPILER OFF CACHE STRING
  "Compiler for matlab/octave interface: mex or mkoctfile or OFF") 
set_property (CACHE MATLAB_COMPILER PROPERTY STRINGS "mex" "mkoctfile" OFF)

# (5) Where to look for data files.  Various classes look in the geoids,
# gravity, magnetic, subdirectories of ${GEOGRAPHICLIB_DATA}.
if (WIN32)
  # The binary installers for the data files for Windows are created
  # with Inno Setup which uses {commonappdata}.  On most Windows
  # systems this is
  # "C:/Documents and Settings/All Users/Application Data", while on
  # newer systems (Windows 7), it is C:/ProgramData.  However the
  # longer name "works" on all Windows systems.
  set (GEOGRAPHICLIB_DATA
    "C:/Documents and Settings/All Users/Application Data/GeographicLib"
    CACHE PATH "Location for data for GeographicLib")
else ()
  set (GEOGRAPHICLIB_DATA
    "/usr/local/share/GeographicLib"
    CACHE PATH "Location for data for GeographicLib")
endif ()

# (6) Compile the short example programs?
if (WIN32)
  option (GEOGRAPHICLIB_EXAMPLES
    "Build example programs in examples directory" OFF)
else ()
  option (GEOGRAPHICLIB_EXAMPLES
    "Build example programs in examples directory" ON)
endif ()

# (7) Create the documentation?  This depends on whether doxygen can be
# found.  If this is OFF, then links will be provided to the online
# documentation on Sourceforge.
if (WIN32)
  option (GEOGRAPHICLIB_DOCUMENTATION
    "Use doxygen to create the documentation" OFF)
else ()
  option (GEOGRAPHICLIB_DOCUMENTATION
    "Use doxygen to create the documentation" ON)
endif ()

# The debug version of the library is called Geographic_d.
set (CMAKE_DEBUG_POSTFIX _d)

# Determine system properties
include (CheckTypeSize)
check_type_size ("long double" LONG_DOUBLE BUILTIN_TYPES_ONLY)
include (TestBigEndian)
test_big_endian (WORDS_BIGENDIAN)

# Create a Config.h to expose system information to the compiler
configure_file (
  include/GeographicLib/Config.h.in
  include/GeographicLib/Config.h )

# The documentation depends on doxygen.  If the documentation is 
if (GEOGRAPHICLIB_DOCUMENTATION)
  set (DOXYGEN_SKIP_DOT ON)
  find_package (Doxygen)
endif ()

# The man pages are written as pod files and converted to nroff format,
# C++ code, and html.  Because this require tools that may not be
# available on an end-user's system, the creation of the final
# documentation is therefore only done in "MAINTAINER" mode.  The
# maintainer runs "make distrib-all" which installs the transformed
# documentation files into the source tree.
if (NOT WIN32)
  find_program (HAVE_POD2MAN pod2man)
  find_program (HAVE_POD2HTML pod2html)
endif ()
if (NOT WIN32 AND HAVE_POD2MAN AND HAVE_POD2HTML)
  set (MAINTAINER ON)
else ()
  set (MAINTAINER OFF)
endif ()
if (MAINTAINER)
  add_custom_target (distrib-all)
  add_dependencies (distrib-all distrib-man)
endif ()

# Look for the tool to compile the Matlab interfaces.
if (MATLAB_COMPILER)
  if (WIN32)
    set (MATLAB_COMPILER_EXT ".bat")
  else ()
    set (MATLAB_COMPILER_EXT "")
  endif ()
  find_program (MEX "${MATLAB_COMPILER}${MATLAB_COMPILER_EXT}")
  if (MATLAB_COMPILER MATCHES "mex")
    get_filename_component (MATLABDIR "${MEX}" REALPATH)
    get_filename_component (MATLABDIR "${MATLABDIR}" PATH)
    find_program (MEXEXTPROG "mexext${MATLAB_COMPILER_EXT}"
      PATHS "${MATLABDIR}")
    execute_process (COMMAND "${MEXEXTPROG}"
      OUTPUT_VARIABLE MEXEXT OUTPUT_STRIP_TRAILING_WHITESPACE)
  else ()
    set (MEXEXT "mex")
  endif ()
  if (NOT MEX)
    message (WARNING
      "Cannot find Matlab compiler ${MATLAB_COMPILER}${MATLAB_COMPILER_EXT}")
  elseif (NOT MEXEXT)
    set (MEX OFF)
    message (WARNING "Cannot determine extension for Matlab compiled code")
  endif ()
endif ()

# Set a default build type for single-configuration cmake generators if
# no build type is set.
if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif ()

# On non-Windows machine, make the compiler more picky.
if (NOT WIN32)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif ()

# The list of tools (to be installed into, e.g., /usr/local/bin)
set (TOOLS CartConvert ConicProj GeodesicProj GeoConvert Geod
  GeoidEval Gravity MagneticField Planimeter TransverseMercatorProj)
# The list of scripts (to be installed into, e.g., /usr/local/sbin)
set (SCRIPTS
  geographiclib-get-geoids geographiclib-get-gravity geographiclib-get-magnetic)

# Set the include directories.  Look in ${PROJECT_BINARY_DIR}/include
# first because that's where Config.h will be
include_directories ("${PROJECT_BINARY_DIR}/include" include)

# The list of subdirectories to process
add_subdirectory (src)
add_subdirectory (include/GeographicLib)
add_subdirectory (tools)
add_subdirectory (man)
add_subdirectory (doc)
add_subdirectory (matlab)
add_subdirectory (python/geographiclib)
if (GEOGRAPHICLIB_EXAMPLES)
  add_subdirectory (examples)
endif ()
add_subdirectory (cmake)

# Packaging support; we deal with 
# (1) a source distribution: cmake make a tar.gz file and the zip file
# is created from this.  Only the maintainer can do this, because of
# the need to generate additional documentation files.
# (2) a binary distribution: code is included for Linux, Apple, and
# Windows, but only the Windows distribution has been exercised.

# Need to ensure that system dlls get included in a binary distribution
if (NOT DEFINED CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_NO_WARNINGS)
  # Visual Studio Express does include redistributable components so
  # squelch the warning.
  set (CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_NO_WARNINGS ON)
endif ()
# To include debug versions of runtime support on Windows, do...
# set (CMAKE_INSTALL_DEBUG_LIBRARIES ON)
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
# The list of files to be excluded from the source distribution.
set (CPACK_SOURCE_IGNORE_FILES
  "#"
  "~\$"
  "/\\\\.git"
  "${PROJECT_SOURCE_DIR}/BUILD"
  "${PROJECT_SOURCE_DIR}/(tests|testdata|cgi-bin|.*\\\\.cache)/"
  "${PROJECT_SOURCE_DIR}/(distrib|.*-distrib|.*-installer|geodesic-papers)/"
  "${PROJECT_SOURCE_DIR}/[^/]*\\\\.(html|kmz)\$"
  "${PROJECT_SOURCE_DIR}/(autogen|biblio|js-compress)\\\\.sh\$"
  "${PROJECT_SOURCE_DIR}/(geodesic-biblio.txt|makefile-admin)\$" )
set (CPACK_SOURCE_GENERATOR TGZ)

set (CPACK_RESOURCE_FILE_LICENSE ${PROJECT_SOURCE_DIR}/LICENSE.txt)
set (CPACK_PACKAGE_INSTALL_DIRECTORY
  "${PROJECT_NAME}-${PACKAGE_VERSION}")
set (CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY}")

if (WIN32)
  # The Windows binary packager is NSIS.  Set the necessary variables
  # for this.
  set (CPACK_NSIS_CONTACT "charles@karney.com")
  set (CPACK_NSIS_URL_INFO_ABOUT "http://geographiclib.sf.net")
  set (CPACK_NSIS_HELP_LINK "mailto:charles@karney.com")
  if (CMAKE_SIZEOF_VOID_P EQUAL 8)
    # Hardcode the prefix for Visual Studio 10
    set (CPACK_NSIS_INSTALL_ROOT "C:\\\\pkg-vc10-x64")
    set (CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY}-win64")
    set (CPACK_NSIS_PACKAGE_NAME
      "${PROJECT_NAME} x64 ${PACKAGE_VERSION}")
    set (CPACK_PACKAGE_INSTALL_REGISTRY_KEY
      "${PROJECT_NAME}-x64-${PACKAGE_VERSION}")
  else ()
    # Hardcode the prefix for Visual Studio 10
    set (CPACK_NSIS_INSTALL_ROOT "C:\\\\pkg-vc10")
    set (CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY}-win32")
    set (CPACK_NSIS_PACKAGE_NAME
      "${PROJECT_NAME} ${PACKAGE_VERSION}")
    set (CPACK_PACKAGE_INSTALL_REGISTRY_KEY
      "${PROJECT_NAME}-${PACKAGE_VERSION}")
  endif ()
  set (CPACK_NSIS_DISPLAY_NAME ${CPACK_NSIS_PACKAGE_NAME})
  set (CPACK_NSIS_MENU_LINKS
    "doc/html/index.html" "Library Documentation"
    "doc/html/utilities.html" "Utilities Documentation"
    "http://geographiclib.sf.net/html/" "Online Documentation")
  set (CPACK_NSIS_MODIFY_PATH ON)
elseif (APPLE)
  # Not tested
  set (CPACK_GENERATOR DMG)
  set (CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY}-darwin")
else ()
  # Not tested
  set (CPACK_GENERATOR TGZ)
endif ()

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
      sed -e "s/Unconfigured/${PACKAGE_VERSION}/"
      include/GeographicLib/Config.h > include/GeographicLib/Config.h.new &&
      mv include/GeographicLib/Config.h.new include/GeographicLib/Config.h
    COMMAND
      cd _CPack_Packages/Linux-Source/TGZ/${CPACK_SOURCE_PACKAGE_FILE_NAME} &&
      touch man/[A-Za-z]*.usage man/[A-Za-z]*.1 man/[A-Za-z]*.1.html &&
      chmod -R g-w .
    COMMAND
      cd _CPack_Packages/Linux-Source/TGZ &&
      find ${CPACK_SOURCE_PACKAGE_FILE_NAME} -type f |
      tar cfzT ${CMAKE_BINARY_DIR}/${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz -
    COMMAND
      rm -f ${CMAKE_BINARY_DIR}/${CPACK_SOURCE_PACKAGE_FILE_NAME}.zip &&
      rsync -a --delete
      _CPack_Packages/Linux-Source/TGZ/${CPACK_SOURCE_PACKAGE_FILE_NAME}
      _CPack_Packages/Linux-Source/TGZ.DOS/ &&
      cd _CPack_Packages/Linux-Source/TGZ.DOS &&
      find . -type f |
      egrep  '/\(doxyfile.in|MANIFEST.in|NEWS|AUTHORS|INSTALL|.*\\.\(cpp|hpp|h.in|txt|pro|usage|pod|py|m|mac|cmake.in|cmake|h|js\)\)$$' |
      xargs unix2dos -q -k &&
      find ${CPACK_SOURCE_PACKAGE_FILE_NAME} -type f |
      zip -q ${CMAKE_BINARY_DIR}/${CPACK_SOURCE_PACKAGE_FILE_NAME}.zip -@
  )
  add_dependencies (dist distrib-all)
endif ()

# Add a test target; the tests are in tools.
enable_testing ()
