# Where the html versions of the man pages (extension .1.html) are
# found.
set (MANDIR ${PROJECT_BINARY_DIR}/man)

# Build up a list of the .1.html files.
set (HTMLMAN)
foreach (TOOL ${TOOLS})
  set (HTMLMAN ${HTMLMAN} ${MANDIR}/${TOOL}.1.html)
endforeach ()

if (COMMON_INSTALL_PATH)
  set (INSTALL_DOC_DIR "share/doc/GeographicLib")
else ()
  set (INSTALL_DOC_DIR "doc")
endif ()

# Run doxygen, if available

# First assemble a list of all the files the documentation uses.  Add a
# dependency on htmlman (from man/CMakeLists.txt).  Use html/index.html
# as the make target.  To make this target, copy the non-doxygen
# generated files into html/.  Run doxfile.in thru cmake's config
# process so that absolute path names are used and so that the pathnames
# are properly stripped by doxygen (via STRIP_FROM_PATH). The
# distrib-doc target copies the html directory into the source tree.

# If doxygen is not available, only the install step (from the source
# tree) is done.

file (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html-stage)
if (DOXYGEN_FOUND)
  set (DOCTARGETS)

  configure_file (GeographicLib.dox.in GeographicLib.dox @ONLY)
  configure_file (doxyfile.in doxyfile @ONLY)
  file (GLOB CXXSOURCES
    ../src/[A-Za-z]*.cpp ../include/GeographicLib/[A-Za-z]*.hpp
    ../tools/[A-Za-z]*.cpp ../examples/[A-Za-z]*.cpp
    ../examples/[A-Za-z]*.hpp)
  file (GLOB EXTRA_FILES ../maxima/[A-Za-z]*.mac
    tmseries30.html geodseries30.html ../LICENSE.txt)
  file (GLOB FIGURES *.png *.svg *.gif)
  file (COPY ${EXTRA_FILES} DESTINATION html-stage)
  add_custom_target (cxxdoc ALL
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/html/index.html)
  add_dependencies (cxxdoc htmlman)
  add_custom_command (OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doxyfile
      ${CMAKE_CURRENT_BINARY_DIR}/GeographicLib.dox
      ${CXXSOURCES} ${EXTRA_FILES} ${FIGURES} ${HTMLMAN}
    COMMAND ${CMAKE_COMMAND} -E remove_directory html
    COMMAND ${CMAKE_COMMAND} -E copy_directory html-stage html
    COMMAND ${DOXYGEN_EXECUTABLE} doxyfile > doxygen.log
    COMMENT "Generating C++ documentation tree")
  set (DOCTARGETS ${DOCTARGETS} cxxdoc)

  add_custom_target (doc ALL)
  add_dependencies (doc ${DOCTARGETS})

  install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html
    DESTINATION ${INSTALL_DOC_DIR})
else ()
  file (COPY ../LICENSE.txt DESTINATION html)
  # Need to absolute path on destination to support old versions of cmake
  configure_file (index.html.in html/index.html)
  configure_file (utilities.html.in html/utilities.html)
  install (FILES
    ${CMAKE_CURRENT_BINARY_DIR}/html/LICENSE.txt
    ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
    ${CMAKE_CURRENT_BINARY_DIR}/html/utilities.html
    DESTINATION ${INSTALL_DOC_DIR}/html)
endif ()
