# $Id$

file (GLOB MATLAB_FILES *.m)
file (GLOB MATLAB_INTERFACE *.cpp)

install (FILES ${MATLAB_FILES} ${MATLAB_INTERFACE}
  DESTINATION share/GeographicLib/matlab)
