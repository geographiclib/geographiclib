# $Id$

# Set up the matlab interface.  The .m files just contain the
# documentation.  The .cpp files implement the interface.
file (GLOB MATLAB_FILES *.m)
file (GLOB MATLAB_INTERFACES *.cpp)
# Both sets of files get installed
install (FILES ${MATLAB_FILES} ${MATLAB_INTERFACES}
  DESTINATION libexec/GeographicLib/matlab)

# If HAVE_MEX then compile the interface routines.  On non-Windows
# systems, an attempt is made to minimize the number of addition target
# generated.  For some reason this does not work with Windows systems.
# On Windows systems, need also to copy the shared library to where the
# mex files live (in the build and install trees).
if (HAVE_MEX)
  if (WIN32)
    include (InstallRequiredSystemLibraries)
    install (PROGRAMS ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS}
      DESTINATION libexec/GeographicLib/matlab)
    add_custom_target (matlab-all)
  else (WIN32)
    set (INTERFACE_LIST)
  endif (WIN32)
  foreach (INTERFACE ${MATLAB_INTERFACES})
    get_filename_component (TARGET ${INTERFACE} NAME_WE)
    set (TARGET_INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}.${MEXEXT})
    if (WIN32)
      add_custom_target (matlab-${TARGET} DEPENDS ${TARGET_INTERFACE})
      add_custom_command (OUTPUT ${TARGET_INTERFACE}
	COMMAND
	  ${MEX}
	    -I${PROJECT_BINARY_DIR}/include -I${PROJECT_SOURCE_DIR}/include
	    -L${PROJECT_BINARY_DIR}/src/Release -lGeographic
	    ${INTERFACE}
	COMMENT "Building matlab interface for ${TARGET}"
	DEPENDS ${INTERFACE} Geographic)
      add_dependencies (matlab-${TARGET} Geographic)
      add_dependencies (matlab-all matlab-${TARGET})
      install (PROGRAMS ${TARGET_INTERFACE}
	DESTINATION libexec/GeographicLib/matlab CONFIGURATIONS Release)
    else (WIN32)
      set (INTERFACE_LIST ${INTERFACE_LIST} ${TARGET_INTERFACE})
      add_custom_command (OUTPUT ${TARGET_INTERFACE}
	COMMAND
	  ${MEX} --mex
	    -I${PROJECT_BINARY_DIR}/include -I${PROJECT_SOURCE_DIR}/include
	    -L${PROJECT_BINARY_DIR}/src -lGeographic
	    -Wl,-rpath=${CMAKE_INSTALL_PREFIX}/lib
	    ${INTERFACE}
	COMMENT "Building matlab interface for ${TARGET}"
	DEPENDS ${INTERFACE} Geographic)
    endif (WIN32)
  endforeach (INTERFACE)
  if (NOT WIN32)
    add_custom_target (matlab-all DEPENDS ${INTERFACE_LIST})
    add_dependencies (matlab-all Geographic)
    install (PROGRAMS ${INTERFACE_LIST}
      DESTINATION libexec/GeographicLib/matlab CONFIGURATIONS Release)
  endif (NOT WIN32)
  if (WIN32 AND GEOGRAPHIC_SHARED_LIB)
    get_target_property (GEOGRAPHIC_LIBRARY_FILE Geographic LOCATION_RELEASE)
    add_custom_command (TARGET matlab-all POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy ${GEOGRAPHIC_LIBRARY_FILE} ./)
    get_filename_component (GEOGRAPHIC_LIBRARY_NAME
      ${GEOGRAPHIC_LIBRARY_FILE} NAME)
    install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${GEOGRAPHIC_LIBRARY_NAME}
      DESTINATION libexec/GeographicLib/matlab CONFIGURATIONS Release)
  endif (WIN32 AND GEOGRAPHIC_SHARED_LIB)
endif (HAVE_MEX)
