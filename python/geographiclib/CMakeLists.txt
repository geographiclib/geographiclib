# $Id$

file (GLOB PYTHON_FILES *.py)
install (FILES ${PYTHON_FILES}
  DESTINATION lib/python/site-packages/geographiclib)
if (WIN32)
  # Only install under Windows because the installation tree is, e.g.,
  # c:/Program Files/GeographicLib/lib/python/site-packages.
  # Don't install on other systems where this command would create,
  # e.g., /usr/local/lib/python/site-packages/setup.py
  install (FILES ../setup.py DESTINATION lib/python/site-packages)
endif (WIN32)
