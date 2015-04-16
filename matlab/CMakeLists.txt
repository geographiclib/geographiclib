# Install matlab files.
if (COMMON_INSTALL_PATH)
  set (INSTALL_MATLAB_DIR "libexec/GeographicLib/matlab")
else ()
  set (INSTALL_MATLAB_DIR "matlab")
endif ()
file (GLOB MATLAB_FILES [A-Za-z]*.m)
install (FILES ${MATLAB_FILES} DESTINATION ${INSTALL_MATLAB_DIR})
# Install "private" functions
file (GLOB PRIVATE_MATLAB_FILES private/[A-Za-z]*.m)
install (FILES ${PRIVATE_MATLAB_FILES}
  DESTINATION ${INSTALL_MATLAB_DIR}/private)
