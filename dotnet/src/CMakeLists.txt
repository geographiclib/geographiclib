# Build the library...

# Include all the .cpp files in the library.
file (GLOB SOURCES [A-Za-z]*.cpp)
file (GLOB HEADERS [A-Za-z]*.h)

add_library (NETGeographic SHARED ${SOURCES} ${HEADERS})
set_target_properties (NETGeographic PROPERTIES COMPILE_FLAGS "/clr")

if (CMAKE_CXX_FLAGS_DEBUG MATCHES "/RTC1")
  string (REPLACE "/RTC1" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
endif ()

if (CMAKE_CXX_FLAGS MATCHES "/EHsc")
  string(REPLACE "/EHsc" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif ()

target_link_libraries (NETGeographic Geographic)

# Set the version number on the library
set_target_properties (NETGeographic PROPERTIES VERSION "${LIBVERSIONFULL}" )

# Specify where the library is installed, adding it to the export depends
install (TARGETS NETGeographic
  EXPORT depends
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib${LIB_SUFFIX}
  ARCHIVE DESTINATION lib${LIB_SUFFIX})

if (MSVC)
  # Install pdb files in debug mode
  get_property(PDB TARGET Geographic PROPERTY LOCATION_DEBUG)
  get_filename_component(_P ${PDB} PATH)
  get_filename_component(PDB ${PDB} NAME_WE)
  set(PDB ${_P}/${PDB}.pdb)
  INSTALL(FILES ${PDB} DESTINATION lib${LIB_SUFFIX} CONFIGURATIONS Debug)
endif()

# Put all the library into a folder in the IDE
set_property (TARGET NETGeographic PROPERTY FOLDER library)
