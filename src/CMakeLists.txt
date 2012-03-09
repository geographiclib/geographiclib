# $Id$

# Build the library...

# Include all the .cpp files in the library.
file (GLOB SOURCES *.cpp)

# Define the library and specify whether it is shared or not.
if (GEOGRAPHIC_SHARED_LIB)
  add_library (Geographic SHARED ${SOURCES})
else (GEOGRAPHIC_SHARED_LIB)
  add_library (Geographic STATIC ${SOURCES})
endif (GEOGRAPHIC_SHARED_LIB)

# Set the version number on the library
if (WIN32)
  set_target_properties (Geographic PROPERTIES VERSION "${LIBVERSIONFULL}" )
else (WIN32)
  set_target_properties (Geographic PROPERTIES
    VERSION "${LIBVERSIONFULL}" SOVERSION "${LIBVERSION}" )
endif (WIN32)

# Specify where the library is installed, adding it to the export DEPENDS
install (TARGETS Geographic
  EXPORT DEPENDS
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
