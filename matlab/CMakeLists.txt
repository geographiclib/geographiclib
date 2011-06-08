# $Id$

file (GLOB MATLAB_FILES *.m)
file (GLOB MATLAB_INTERFACES *.cpp)
install (FILES ${MATLAB_FILES} ${MATLAB_INTERFACES}
  DESTINATION share/GeographicLib/matlab)

if (HAVE_MEX)
  add_custom_target (MATLAB-all ALL)
  foreach (INTERFACE ${MATLAB_INTERFACES})
    get_filename_component (TARGET ${INTERFACE} NAME_WE)
    set (TARGET_INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}.mexw32)
    add_custom_target (MATLAB-${TARGET} DEPENDS ${TARGET_INTERFACE})
    add_custom_command (OUTPUT ${TARGET_INTERFACE}
      COMMAND
	${HAVE_MEX}
	  -I${PROJECT_BINARY_DIR}/include -I${PROJECT_SOURCE_DIR}/include
	  -L${PROJECT_BINARY_DIR}/src/Release -lGeographic 
	  ${INTERFACE}
      DEPENDS ${INTERFACE} Geographic)
    add_dependencies (MATLAB-all MATLAB-${TARGET})
    install (FILES ${TARGET_INTERFACE}
      DESTINATION share/GeographicLib/matlab
      CONFIGURATIONS Release)
  endforeach (INTERFACE)
endif (HAVE_MEX)
