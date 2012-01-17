# $Id$

# Where the html versions of the man pages (extension .1.html) are
# found.
if (MAINTAINER)
  set (MANDIR ${PROJECT_BINARY_DIR}/man)
else (MAINTAINER)
  set (MANDIR ${PROJECT_SOURCE_DIR}/man)
endif (MAINTAINER)

# BUILD UP A LIST OF THE .1.HTML FILES.
set (HTMLMAN)
foreach (TOOL ${TOOLS})
  set (HTMLMAN ${HTMLMAN} ${MANDIR}/${TOOL}.1.html)
endforeach (TOOL)

# The maintainer runs doxygen...

# First assemble a list of all the files the documentation uses.  Add a
# dependency on htmlman (from man/CMakeLists.txt).  Use doxygen.log as
# the make target.  To make this target, copy the non-doxygen generated
# files into html/.  Run Doxfile thru a sed command so that absolute
# path names are used and so that the pathnames are properly stripped by
# doxygen (via STRIP_FROM_PATH). The distrib-doc target copies the html
# directory into the source tree.

# In non-maintenance mode, only the install step (from the source tree)
# is needed.

if (MAINTAINER)
  file (GLOB SOURCES
    ../src/*.cpp ../include/GeographicLib/*.hpp
    ../tools/*.cpp ../examples/*.cpp)
  file (GLOB EXTRA_FILES ../maxima/*.mac tmseries30.html geodseries30.html
    ../LICENSE.txt)
  add_custom_target (doc ALL
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doxygen.log)
  add_dependencies (doc htmlman)
  add_custom_command (OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/doxygen.log
    DEPENDS Doxyfile Geographic.doc ${SOURCES} ${EXTRA_FILES} ${HTMLMAN}
    COMMAND rm -rf html && mkdir html
    COMMAND
      cp -p ${EXTRA_FILES} ${HTMLMAN} html/
    COMMAND
      sed \"s%\\.\\./%${PROJECT_SOURCE_DIR}/%g\"
        ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile | doxygen - > doxygen.log
    COMMENT "Generating html documentation tree")
  install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html
    DESTINATION share/doc/GeographicLib)

  add_custom_target (distrib-doc)
  add_dependencies (distrib-doc doc)
  add_custom_command (TARGET distrib-doc
    COMMAND rsync -a --delete html/ ${CMAKE_CURRENT_SOURCE_DIR}/html/
    COMMAND touch ${CMAKE_CURRENT_SOURCE_DIR}/html/index.html
    COMMENT "Installing html documentation tree in source tree")

else (MAINTAINER)
  install (DIRECTORY html DESTINATION share/doc/GeographicLib)
endif (MAINTAINER)

# Finally install the Javascript files
file (GLOB SCRIPTDRIVERS scripts/*.html)
file (GLOB JSSCRIPTS scripts/GeographicLib/*.js)
install (FILES ${SCRIPTDRIVERS} DESTINATION share/doc/GeographicLib/scripts)
install (FILES ${JSSCRIPTS}
  DESTINATION share/doc/GeographicLib/scripts/GeographicLib)
