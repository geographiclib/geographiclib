# Build the .NET library...

# Include all the .cpp files in the library.
file (GLOB SOURCES [A-Za-z]*.cpp)
file (GLOB HEADERS [A-Za-z]*.h)

add_library (${NETGEOGRAPHICLIB_LIBRARIES} SHARED ${SOURCES} ${HEADERS})
set_target_properties (${NETGEOGRAPHICLIB_LIBRARIES}
  PROPERTIES COMPILE_FLAGS "/clr")

string (REPLACE "/RTC1" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
string (REPLACE "/EHsc" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

add_definitions (${PROJECT_DEFINITIONS})
target_link_libraries (${NETGEOGRAPHICLIB_LIBRARIES} ${PROJECT_LIBRARIES})

# Set the version number on the library
set_target_properties (${NETGEOGRAPHICLIB_LIBRARIES} PROPERTIES
  VERSION "${LIBVERSIONFULL}" OUTPUT_NAME ${NETLIBNAME})

# Specify where the library is installed, adding it to the export depends
install (TARGETS ${NETGEOGRAPHICLIB_LIBRARIES}
  EXPORT depends
  RUNTIME DESTINATION bin)

if (PACKAGE_DEBUG_LIBS)
  install (PROGRAMS
    "${PROJECT_BINARY_DIR}/bin/Debug/${NETLIBNAME}${CMAKE_DEBUG_POSTFIX}.dll"
    DESTINATION bin CONFIGURATIONS Release)
endif ()

# Install pdb file.
foreach (_c ${CMAKE_CONFIGURATION_TYPES})
  string (TOUPPER ${_c} _C)
  get_target_property (_P ${NETGEOGRAPHICLIB_LIBRARIES} LOCATION_${_C})
  get_filename_component (_D ${_P} PATH)
  get_filename_component (_N ${_P} NAME_WE)
  set (_PDB ${_D}/${_N}.pdb)
  install (FILES ${_PDB} DESTINATION bin CONFIGURATIONS ${_c} OPTIONAL)
endforeach ()

# Put all the library into a folder in the IDE
set_property (TARGET ${NETGEOGRAPHICLIB_LIBRARIES} PROPERTY FOLDER library)
