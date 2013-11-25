# Build the .NET library...

# By default this is only invoked for Visual Studio generators other
# than Visual Studio 10.  Visual Studio 10 generates triggers a bunch of
# warnings of the form:
#
# 3>  Geodesic.cpp
# 3>c:\program files (x86)\reference assemblies\microsoft\framework\.netframework\v4.0\system.core.dll : warning C4945: 'ExtensionAttribute' : cannot import symbol from 'c:\program files (x86)\reference assemblies\microsoft\framework\.netframework\v4.0\system.core.dll': as 'System::Runtime::CompilerServices::ExtensionAttribute' has already been imported from another assembly 'mscorlib'
# 3> c:\windows\microsoft.net\framework\v4.0.30319\mscorlib.dll : see declaration of 'System::Runtime::CompilerServices::ExtensionAttribute'
# 3> first seen type is used; re-order imported assemblies to use the current type
# 3> This diagnostic occurred while importing type 'System.Runtime.CompilerServices.ExtensionAttribute' from assembly 'System.Core, Version=4.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089'.
#
# It appears however that these warnings can be safely ignored.

# Include all the .cpp files in the library.
file (GLOB SOURCES [A-Za-z]*.cpp)
file (GLOB HEADERS [A-Za-z]*.h)

add_library (${NETGEOGRAPHICLIB_LIBRARIES} SHARED ${SOURCES} ${HEADERS})
set_target_properties (${NETGEOGRAPHICLIB_LIBRARIES}
  PROPERTIES COMPILE_FLAGS "/clr")

string (REPLACE "/RTC1" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
string (REPLACE "/EHsc" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

target_link_libraries (${NETGEOGRAPHICLIB_LIBRARIES} ${PROJECT_LIBRARIES})

# Set the version number on the library
set_target_properties (${NETGEOGRAPHICLIB_LIBRARIES} PROPERTIES
  VERSION "${LIBVERSIONFULL}" PDB_NAME NETGeographicLib)

# Specify where the library is installed, adding it to the export depends
install (TARGETS ${NETGEOGRAPHICLIB_LIBRARIES}
  EXPORT depends
  RUNTIME DESTINATION bin)

if (MSVC AND PACKAGE_DEBUG_LIBS)
  install (PROGRAMS
    "${CMAKE_CURRENT_BINARY_DIR}/Debug/NETGeographicLib${CMAKE_DEBUG_POSTFIX}.dll"
    DESTINATION bin CONFIGURATIONS Release)
endif ()

if (MSVC AND NOT CMAKE_VERSION VERSION_LESS 2.8.11)
  # Install pdb files in debug mode.  Names of PDB files are wrong for
  # cmake 2.8.0 (presumably because the PDB_NAME property on the library
  # is ignored; so require version >= 2.8.11 (but this is likely to be
  # too conservative).
  get_target_property (_P NETGeographicLib LOCATION_DEBUG)
  get_filename_component (_P ${_P} PATH)
  get_target_property (_N NETGeographicLib PDB_NAME)
  set (_PDB ${_P}/${_N}.pdb)
  install (FILES ${_PDB} DESTINATION bin
    RENAME NETGeographicLib${CMAKE_DEBUG_POSTFIX}.pdb CONFIGURATIONS Debug)
endif()

# Put all the library into a folder in the IDE
set_property (TARGET ${NETGEOGRAPHICLIB_LIBRARIES} PROPERTY FOLDER library)
