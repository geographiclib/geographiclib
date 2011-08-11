# $Id$

file (GLOB MATLAB_FILES *.m)
file (GLOB MATLAB_INTERFACES *.cpp)
install (FILES ${MATLAB_FILES} ${MATLAB_INTERFACES}
  DESTINATION libexec/GeographicLib/matlab)

if (HAVE_MEX)
  if (WIN32)
    include (InstallRequiredSystemLibraries)
    install (PROGRAMS ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS}
      DESTINATION libexec/GeographicLib/matlab)
  endif (WIN32)
  add_custom_target (matlab-all ALL)
  foreach (INTERFACE ${MATLAB_INTERFACES})
    get_filename_component (TARGET ${INTERFACE} NAME_WE)
    set (TARGET_INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}.${MEXEXT})
    add_custom_target (matlab-${TARGET} DEPENDS ${TARGET_INTERFACE})
    if (WIN32)
      add_custom_command (OUTPUT ${TARGET_INTERFACE}
	COMMAND
	  ${MEX}
	    -I${PROJECT_BINARY_DIR}/include -I${PROJECT_SOURCE_DIR}/include
	    -L${PROJECT_BINARY_DIR}/src/${CMAKE_CFG_INTDIR} -lGeographic
	    ${INTERFACE}
	COMMENT "Building matlab interface for ${TARGET}"
	DEPENDS ${INTERFACE} Geographic)
    else (WIN32)
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
    add_dependencies (matlab-all matlab-${TARGET})
    install (PROGRAMS ${TARGET_INTERFACE}
      DESTINATION libexec/GeographicLib/matlab)
  endforeach (INTERFACE)
  if (WIN32 AND GEOGRAPHIC_SHARED_LIB)
    get_target_property (GEOGRAPHIC_LIBRARY_FILE Geographic
      "LOCATION_${CMAKE_CFG_INTDIR}")
    add_custom_command (TARGET matlab-all POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy ${GEOGRAPHIC_LIBRARY_FILE} ./)
    get_filename_component (GEOGRAPHIC_LIBRARY_NAME
      ${GEOGRAPHIC_LIBRARY_FILE} NAME)
    install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${GEOGRAPHIC_LIBRARY_NAME}
      DESTINATION libexec/GeographicLib/matlab)
  endif (WIN32 AND GEOGRAPHIC_SHARED_LIB)
endif (HAVE_MEX)
