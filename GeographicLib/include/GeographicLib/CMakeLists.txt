# $Id$

file (GLOB HEADERS *.hpp)
install (FILES ${HEADERS} DESTINATION include/GeographicLib)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/Config.h DESTINATION include/GeographicLib)
