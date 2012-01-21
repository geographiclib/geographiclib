# $Id: ad1a97beeda9134b93debee2ba855548e252c3a7 $

# Install the header files including the generated Config.h (which is in
# the build tree).
file (GLOB HEADERS *.hpp)
install (FILES ${HEADERS} DESTINATION include/GeographicLib)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/Config.h
  DESTINATION include/GeographicLib)
