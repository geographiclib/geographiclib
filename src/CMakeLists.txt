# $Id$

file (GLOB SOURCES *.cpp)

if (GEOGRAPHIC_SHARED_LIB)
  add_library (Geographic SHARED ${SOURCES})
else (GEOGRAPHIC_SHARED_LIB)
  add_library (Geographic STATIC ${SOURCES})
endif (GEOGRAPHIC_SHARED_LIB)

target_link_libraries (Geographic)
if (WIN32)
  set_target_properties (Geographic PROPERTIES VERSION "${LIBVERSIONFULL}" )
else (WIN32)
  set_target_properties (Geographic PROPERTIES
    VERSION "${LIBVERSIONFULL}" SOVERSION "${LIBVERSION}" )
endif (WIN32)
install (TARGETS Geographic
  EXPORT GeographicLibDepends
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
