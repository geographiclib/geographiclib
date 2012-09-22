# config file support for find_package(GeographicLib).  This needs to
# deal with two environments: (1) finding the build tree and (2)
# finding the install tree.  geographiclib-config.cmake detects which
# situation it is handing by looking at @GEOGRAPHICLIB_ROOT_DIR@.  If
# this is an absolute path, it's in the build tree; otherwise, it's in the
# install tree.  (Note that the whole install tree can be relocated.)

# geographiclib-config.cmake for the build tree
set (GEOGRAPHICLIB_ROOT_DIR "${PROJECT_BINARY_DIR}")
set (GeographicLib_INCLUDE_DIRS
  "${PROJECT_BINARY_DIR}/include" "${PROJECT_SOURCE_DIR}/include" )
configure_file (geographiclib-config.cmake.in
  "${PROJECT_BINARY_DIR}/geographiclib-config.cmake" @ONLY)
configure_file (geographiclib-config-version.cmake.in
  "${PROJECT_BINARY_DIR}/geographiclib-config-version.cmake" @ONLY)
export (TARGETS Geographic ${TOOLS}
  FILE "${PROJECT_BINARY_DIR}/geographiclib-depends.cmake")

# geographiclib-config.cmake for the install tree.  It's installed in
# ${INSTALL_CMAKE_DIR} and @GEOGRAPHICLIB_ROOT_DIR@ is the relative
# path to the root from there.  (Note that the whole install tree can
# be relocated.)
if (COMMON_INSTALL_PATH)
  set (INSTALL_CMAKE_DIR "share/cmake/${PROJECT_NAME}-${PACKAGE_VERSION}")
  set (GEOGRAPHICLIB_ROOT_DIR "../../..")
else ()
  set (INSTALL_CMAKE_DIR "cmake")
  set (GEOGRAPHICLIB_ROOT_DIR "..")
endif ()
# @GeographicLib_INCLUDE_DIRS@ is not used in the install tree; reset
# it to prevent the source and build paths appearing in the installed
# config files
set (GeographicLib_INCLUDE_DIRS IGNORE)
# Need to absolute path on destination to support old versions of cmake
configure_file (geographiclib-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/geographiclib-config.cmake @ONLY)
configure_file (geographiclib-config-version.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/geographiclib-config-version.cmake @ONLY)
install (FILES
  "${CMAKE_CURRENT_BINARY_DIR}/geographiclib-config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/geographiclib-config-version.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}")
# Make information about the cmake targets (the library and the tools)
# available.
install (EXPORT depends
  FILE geographiclib-depends.cmake
  DESTINATION "${INSTALL_CMAKE_DIR}")
