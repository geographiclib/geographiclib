# Build the library...

# Include all the .cpp files in the library.
file (GLOB SOURCES [A-Za-z]*.cpp)
file (GLOB HEADERS
  ${PROJECT_BINARY_DIR}/include/GeographicLib/Config.h
  ../include/GeographicLib/[A-Za-z]*.hpp)

# Define the library and specify whether it is shared or not.
if (GEOGRAPHIC_SHARED_LIB)
  add_library (Geographic SHARED ${SOURCES} ${HEADERS})
else ()
  add_library (Geographic STATIC ${SOURCES} ${HEADERS})
endif ()

# Set the version number on the library
if (WIN32)
  set_target_properties (Geographic PROPERTIES VERSION "${LIBVERSIONFULL}" )
else ()
  set_target_properties (Geographic PROPERTIES
    VERSION "${LIBVERSIONFULL}" SOVERSION "${LIBVERSION}" )
endif ()

# Specify where the library is installed, adding it to the export depends
install (TARGETS Geographic
  EXPORT depends
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

# Put all the library into a folder in the IDE
set_property (TARGET Geographic PROPERTY FOLDER library)
