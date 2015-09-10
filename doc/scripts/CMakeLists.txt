# This list governs the order in which the JavaScript sources are
# concatenated.  This shouldn't be changed.
set (JS_MODULES Math Geodesic GeodesicLine PolygonArea DMS)

# Combine JavaScript into a single file if necessary
set (JS_BUILD 0)
set (JSSCRIPTS)
set (JS_TARGET "${CMAKE_CURRENT_BINARY_DIR}/geographiclib.js")
set (FILE_INVENTORY "")
foreach (_F ${JS_MODULES})
  set (_S "GeographicLib/${_F}.js")
  set (FILE_INVENTORY "${FILE_INVENTORY} ${_F}.js")
  list (APPEND JSSCRIPTS  ${_S})
  if ("${CMAKE_CURRENT_SOURCE_DIR}/_S" IS_NEWER_THAN ${JS_TARGET})
    set (JS_BUILD 1)
  endif ()
endforeach ()

if (JS_BUILD)
  file (STRINGS "GeographicLib/Math.js" _S REGEX version_string)
  string (REGEX REPLACE ".*\"(.*)\".*" "\\1" JS_VERSION "${_S}")
  file (REMOVE ${JS_TARGET})
  file (READ "HEADER.js" _S)
  string (CONFIGURE ${_S} _S @ONLY)
  file (APPEND ${JS_TARGET} "${_S}")
  foreach (_F ${JSSCRIPTS})
    file (READ "${_F}" _S)
    # Normalize the line endings.
    string (REGEX REPLACE "\r" "\n" _S "${_S}")
    # This matches /*...*/ style comments, where ... is any number of
    # \*[^/] and [^*].  This has the defect that the it won't detect,
    # e.g., **/ as the end of the comment.
    string (REGEX REPLACE "/\\*(\\*[^/]|[^*])*\\*/" "" _S "${_S}")
    string (REGEX REPLACE "//[^\n]*\n" "\n" _S "${_S}")
    string (REGEX REPLACE "[ \t]+" " " _S "${_S}")
    string (REGEX REPLACE "([^\"A-Za-z0-9_]) " "\\1" _S "${_S}")
    string (REGEX REPLACE " ([^\\[\"A-Za-z0-9_])"  "\\1" _S "${_S}")
    string (REGEX REPLACE "\n " "\n" _S "${_S}")
    string (REGEX REPLACE " \n" "\n" _S "${_S}")
    string (REGEX REPLACE "^\n" "" _S "${_S}")
    string (REGEX REPLACE "\n+" "\n" _S "${_S}")
    file (APPEND ${JS_TARGET} "// ${_F}\n${_S}")
  endforeach ()
  # support loading with node's require
  file (APPEND ${JS_TARGET}
    "if(typeof module!=='undefined')module.exports=GeographicLib;\n")
endif ()

# "make javascript" will reconfigure cmake if necessary, since
# geographiclib.js is created during configuration.
add_custom_command (OUTPUT geographiclib.js
  DEPENDS ${JSSCRIPTS} HEADER.js
  COMMAND ${CMAKE_COMMAND} ARGS "."
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
add_custom_target (javascript ALL DEPENDS geographiclib.js)

if (COMMON_INSTALL_PATH)
  set (INSTALL_SCRIPTS_DIR "share/doc/GeographicLib/scripts")
else ()
  set (INSTALL_SCRIPTS_DIR "doc/scripts")
endif ()

# Install the JavaScript file
install (FILES ${JS_TARGET} DESTINATION ${INSTALL_SCRIPTS_DIR})

set (SCRIPTDRIVERS geod-calc geod-google geod-google-instructions)
foreach (_D ${SCRIPTDRIVERS})
  configure_file (${_D}.in ${_D}.html)
  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/${_D}.html
    DESTINATION ${INSTALL_SCRIPTS_DIR})
endforeach ()
