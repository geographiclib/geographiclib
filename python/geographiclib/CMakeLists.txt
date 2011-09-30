# $Id: e5c08037ac4e61ea6e2f0c05425deb4495f01bce $

file (GLOB PYTHON_FILES *.py)
install (FILES ${PYTHON_FILES}
  DESTINATION lib/python/site-packages/geographiclib)
install (FILES ../setup.py DESTINATION lib/python/site-packages)
