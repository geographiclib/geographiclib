# Build the tools...

# Where to find the *.usage files for the --help option.
if (RELEASE)
  set (MANDIR ${PROJECT_SOURCE_DIR}/man)
else ()
  set (MANDIR ${PROJECT_BINARY_DIR}/man)
endif ()
include_directories (${MANDIR})

# Only needed if target_compile_definitions is not supported
add_definitions (${PROJECT_DEFINITIONS})

# Loop over all the tools, specifying the source and library.
if (BINDIR)
  add_custom_target (${GEOGRAPHICLIB_TOOLS_TARGET} ALL)
else ()
  add_custom_target (${GEOGRAPHICLIB_TOOLS_TARGET})
endif ()

foreach (TOOL ${TOOLS})

  if (BINDIR)
    add_executable (${TOOL} ${TOOL}.cpp)
  else ()
    add_executable (${TOOL} EXCLUDE_FROM_ALL ${TOOL}.cpp)
  endif ()
  if (NOT RELEASE)
    add_dependencies (${TOOL} usage)
  endif ()
  add_dependencies (${GEOGRAPHICLIB_TOOLS_TARGET} ${TOOL})

  set_source_files_properties (${TOOL}.cpp PROPERTIES
    OBJECT_DEPENDS ${MANDIR}/${TOOL}.usage)

  target_link_libraries (${TOOL} ${PROJECT_LIBRARIES} ${HIGHPREC_LIBRARIES})

endforeach ()

if (MSVC OR CMAKE_CONFIGURATION_TYPES)
  # Add _d suffix for your debug versions of the tools
  set_target_properties (${TOOLS} PROPERTIES
    DEBUG_POSTFIX "${CMAKE_DEBUG_POSTFIX}")
endif ()

if (APPLE AND RELATIVE_LIBDIR)
  # Ensure that the package is relocatable
  set_target_properties (${TOOLS} PROPERTIES
    INSTALL_RPATH "@loader_path/${RELATIVE_LIBDIR}")
endif ()

# Specify where the tools are installed, adding them to the export targets
if (BINDIR)
  install (TARGETS ${TOOLS} EXPORT targets DESTINATION ${BINDIR})
endif ()

# Put all the tools into a folder in the IDE
set_property (TARGET ${GEOGRAPHICLIB_TOOLS_TARGET} ${TOOLS}
  PROPERTY FOLDER tools)

# Create the scripts for downloading the data files on non-Windows
# systems.  This needs to substitute ${GEOGRAPHICLIB_DATA} as the
# default data directory.  These are installed under sbin, because it is
# expected to be run with write access to /usr/local.
if (NOT CMAKE_HOST_WIN32)
  foreach (SCRIPT ${SCRIPTS})
    configure_file (${SCRIPT}.sh scripts/${SCRIPT} @ONLY)
    add_custom_command (OUTPUT ${SCRIPT}
      COMMAND ${CMAKE_COMMAND} -E
        copy scripts/${SCRIPT} ${SCRIPT} && chmod +x ${SCRIPT}
      DEPENDS ${SCRIPT}.sh)
    if (SBINDIR)
      install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}
        DESTINATION ${SBINDIR})
    endif ()
  endforeach ()
  add_custom_target (${GEOGRAPHICLIB_SCRIPTS_TARGET} ALL DEPENDS ${SCRIPTS})
endif ()
