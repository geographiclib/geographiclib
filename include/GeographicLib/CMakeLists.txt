# $Id: 43e619c4d3b8b2ce7857575bfa2029adb739f68e $

file (GLOB HEADERS *.hpp)
install (FILES ${HEADERS} DESTINATION include/GeographicLib)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/Config.h DESTINATION include/GeographicLib)
