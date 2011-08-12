# $Id$

file (GLOB PYTHON_FILES *.py)
install (FILES ${PYTHON_FILES}
  DESTINATION lib/python/site-packages/geographiclib)
