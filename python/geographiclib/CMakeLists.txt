# $Id: 4693af55b838c9d355c9233558af2a71eed8ff00 $

file (GLOB PYTHON_FILES *.py)
install (FILES ${PYTHON_FILES}
  DESTINATION lib/python/site-packages/geographiclib)
