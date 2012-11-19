# Set up the matlab interface.  The .m files just contain the
# documentation.  The .cpp files implement the interface.
file (GLOB MATLAB_FILES [A-Za-z]*.m)
file (GLOB MATLAB_INTERFACES [A-Za-z]*.cpp)
if (COMMON_INSTALL_PATH)
  set (INSTALL_MATLAB_DIR "libexec/GeographicLib/matlab")
else ()
  set (INSTALL_MATLAB_DIR "matlab")
endif ()
# Both sets of files get installed
install (FILES ${MATLAB_FILES} ${MATLAB_INTERFACES}
  DESTINATION ${INSTALL_MATLAB_DIR})
# Install "private" functions
file (GLOB PRIVATE_MATLAB_FILES private/[A-Za-z]*.m)
install (FILES ${PRIVATE_MATLAB_FILES}
  DESTINATION ${INSTALL_MATLAB_DIR}/private)

# If MEX then compile the interface routines.  On non-Windows
# systems, an attempt is made to minimize the number of addition target
# generated.  For some reason this does not work with Windows systems.
# On Windows systems, need also to copy the shared library to where the
# mex files live (in the build and install trees).
if (MEX)
  if (WIN32)
    include (InstallRequiredSystemLibraries)
    install (PROGRAMS ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS}
      DESTINATION ${INSTALL_MATLAB_DIR})
    add_custom_target (matlab-all)
  else ()
    set (INTERFACE_LIST)
  endif ()
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
        DESTINATION ${INSTALL_MATLAB_DIR} CONFIGURATIONS Release)
      # Put all the matlab targets into a folder in the IDE
      set_property (TARGET matlab-${TARGET} PROPERTY FOLDER matlab)
    else ()
      set (INTERFACE_LIST ${INTERFACE_LIST} ${TARGET_INTERFACE})
      add_custom_command (OUTPUT ${TARGET_INTERFACE}
        COMMAND
          ${MEX} --mex
            -I${PROJECT_BINARY_DIR}/include -I${PROJECT_SOURCE_DIR}/include
            -L${PROJECT_BINARY_DIR}/src -lGeographic
            -Wl,-rpath=${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}
            ${INTERFACE}
        COMMENT "Building matlab interface for ${TARGET}"
        DEPENDS ${INTERFACE} Geographic)
    endif ()
  endforeach ()
  if (NOT WIN32)
    add_custom_target (matlab-all DEPENDS ${INTERFACE_LIST})
    add_dependencies (matlab-all Geographic)
    install (PROGRAMS ${INTERFACE_LIST}
      DESTINATION ${INSTALL_MATLAB_DIR} CONFIGURATIONS Release)
  endif ()
  if (WIN32 AND GEOGRAPHIC_SHARED_LIB)
    get_target_property (GEOGRAPHIC_LIBRARY_FILE Geographic LOCATION_RELEASE)
    add_custom_command (TARGET matlab-all POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy ${GEOGRAPHIC_LIBRARY_FILE} ./)
    get_filename_component (GEOGRAPHIC_LIBRARY_NAME
      ${GEOGRAPHIC_LIBRARY_FILE} NAME)
    install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${GEOGRAPHIC_LIBRARY_NAME}
      DESTINATION ${INSTALL_MATLAB_DIR} CONFIGURATIONS Release)
  endif ()

  set_property (TARGET matlab-all PROPERTY FOLDER matlab)
endif ()
