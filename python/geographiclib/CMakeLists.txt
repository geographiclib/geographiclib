# $Id$

# Install the python files.   Probably full-time python users should get
# the python package from http://pypi.python.org/pypi/geographiclib and
# "install" this.
file (GLOB PYTHON_FILES *.py)
install (FILES ${PYTHON_FILES}
  DESTINATION lib/python/site-packages/geographiclib)

# Only install setup.py under Windows because the installation tree is,
# e.g., c:/Program Files/GeographicLib/lib/python/site-packages.  Don't
# install on other systems where this command would create, e.g.,
# /usr/local/lib/python/site-packages/setup.py which is likely to cause
# conflicts.
if (WIN32)
  install (FILES ../setup.py DESTINATION lib/python/site-packages)
endif (WIN32)
