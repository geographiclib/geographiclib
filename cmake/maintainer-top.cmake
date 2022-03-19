set (DISTRIB_DIR "${CMAKE_BINARY_DIR}/distrib")
set (DISTRIB_NAME "${DISTRIB_DIR}/${PACKAGE_DIR}")
add_custom_target (prep-source
  COMMAND ${CMAKE_MAKE_PROGRAM} package_source
  COMMAND ${CMAKE_COMMAND} -E rm -rf ${DISTRIB_NAME}
  COMMAND ${CMAKE_COMMAND} -E copy_directory
  _CPack_Packages/Linux-Source/TGZ/${PACKAGE_DIR} ${DISTRIB_NAME}
  COMMAND cd ${DISTRIB_NAME} &&
  find * -type f | sort -u > ${DISTRIB_DIR}/files.1 &&
  ( cd ${PROJECT_SOURCE_DIR} && git ls-tree -r HEAD --name-only ) |
  sort -u > ${DISTRIB_DIR}/files.2 &&
  comm -23 ${DISTRIB_DIR}/files.[12] | xargs -r -d '\\n' rm
  COMMAND ${CMAKE_COMMAND} -E rm -f autogen.done)
add_custom_command (OUTPUT autogen.done
  COMMAND cd ${DISTRIB_NAME} && ${PROJECT_SOURCE_DIR}/autogen.sh &&
  touch ${PROJECT_BINARY_DIR}/autogen.done
  DEPENDS prep-source autogen.sh configure.ac
  Makefile.am src/Makefile.am include/Makefile.am tools/Makefile.am
  doc/Makefile.am man/Makefile.am cmake/Makefile.am
  examples/Makefile.am tests/Makefile.am)
add_dependencies (distrib-man prep-source)
add_custom_target (distrib-all DEPENDS distrib-man autogen.done)
add_custom_command (TARGET distrib-all
  COMMAND cd ${DISTRIB_NAME} && echo ${PROJECT_VERSION} > VERSION &&
  chmod -R g-w .)
add_custom_target (dist
  COMMAND
  cd ${DISTRIB_DIR} &&
  find ${PACKAGE_DIR} -type f | tar cfzT ${PACKAGE_NAME}.tar.gz -
  COMMAND
  rm -f ${DISTRIB_DIR}/${PACKAGE_NAME}.zip &&
  cd ${DISTRIB_DIR} &&
  find ${PACKAGE_DIR} -type f | zip -q ${PACKAGE_NAME}.zip -@
  COMMENT
  "created distrib/${PACKAGE_NAME}.{tar.gz,zip}")
add_dependencies (dist distrib-all)

if (RSYNC)
  set (USER karney)
  set (DATAROOT $ENV{HOME}/web/geographiclib-files/distrib-C++)
  set (DOCROOT $ENV{HOME}/web/geographiclib-web/htdocs)
  add_custom_target (stage-dist
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${DISTRIB_DIR}/${PACKAGE_NAME}.tar.gz
    ${DISTRIB_DIR}/${PACKAGE_NAME}.zip
    ${PROJECT_SOURCE_DIR}/data-distrib/distrib-C++/
    COMMAND ${RSYNC} --delete -av
    ${PROJECT_SOURCE_DIR}/data-distrib/distrib-C++/ ${DATAROOT}/)
  add_dependencies(stage-dist dist)

  if (BUILD_DOCUMENTATION AND DOXYGEN_FOUND)
    add_custom_target (stage-doc
      COMMAND ${RSYNC} --delete -av
      doc/html/ ${DOCROOT}/C++/${PROJECT_VERSION}${PROJECT_VERSION_SUFFIX}/)
    add_dependencies(stage-doc doc)
  endif ()

  add_custom_target (deploy-dist
    COMMAND ${RSYNC} --delete -av ${DATAROOT}
    ${USER}@frs.sourceforge.net:/home/frs/project/geographiclib/)
  add_custom_target (deploy-doc
    COMMAND ${RSYNC} --delete -av -e ssh
    ${DOCROOT}/C++ ${USER},geographiclib@web.sourceforge.net:./htdocs/)
endif ()

if (NOT WIN32)
  set (BINARY_EXT "m4|gif|pdf|png|kmz")
  add_custom_target (checktrailingspace
    COMMAND git ls-tree -r HEAD --name-only |
    egrep -v '\\.\(${BINARY_EXT}\)$$' |
    xargs grep '[ \t]$$' || true
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Looking for trailing spaces")
  add_custom_target (checktabs
    COMMAND git ls-tree -r HEAD --name-only |
    egrep -v '\([Mm]akefile|\\.\(${BINARY_EXT}\)$$\)' |
    xargs grep -l '\t' || true
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Looking for tabs")
  add_custom_target (checkblanklines
    COMMAND git ls-tree -r HEAD --name-only |
    egrep -v '\\.\(${BINARY_EXT}\)$$' |
    while read f\; do tr 'X\\n' 'YX' < $$f |
    egrep '\(^X|XXX|XX$$|[^X]$$\)' > /dev/null && echo $$f\; done || true
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Looking for extra blank lines")

  add_custom_target(sanitize)
  add_dependencies(sanitize checktrailingspace checktabs checkblanklines)
endif ()
