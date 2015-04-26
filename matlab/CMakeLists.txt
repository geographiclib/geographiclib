# Install matlab files.
if (COMMON_INSTALL_PATH)
  # More Octave friendly would be "share/octave/site/m"
  set (INSTALL_MATLAB_DIR "share/matlab")
else ()
  set (INSTALL_MATLAB_DIR "matlab")
endif ()
file (GLOB MATLAB_FILES geographiclib/[A-Za-z]*.m)
install (FILES ${MATLAB_FILES} DESTINATION ${INSTALL_MATLAB_DIR}/geographiclib)
# Install "private" functions
file (GLOB PRIVATE_MATLAB_FILES geographiclib/private/[A-Za-z]*.m)
install (FILES ${PRIVATE_MATLAB_FILES}
  DESTINATION ${INSTALL_MATLAB_DIR}/geographiclib/private)
# Install "legacy" functions
file (GLOB LEGACY_FILES
  geographiclib_legacy/[A-Za-z]*.m  geographiclib_legacy/[A-Za-z]*.cpp)
install (FILES ${LEGACY_FILES}
  DESTINATION ${INSTALL_MATLAB_DIR}/geographiclib_legacy)
# Octave doesn't know about Contents.m, so
install (FILES geographiclib/Contents.m
  DESTINATION ${INSTALL_MATLAB_DIR}/geographiclib RENAME geographiclib.m)
install (FILES geographiclib_legacy/Contents.m
  DESTINATION ${INSTALL_MATLAB_DIR}/geographiclib_legacy
  RENAME geographiclib_legacy.m)
