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

add_library (NETGeographic SHARED ${SOURCES} ${HEADERS})
set_target_properties (NETGeographic PROPERTIES COMPILE_FLAGS "/clr")

string (REPLACE "/RTC1" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
string (REPLACE "/EHsc" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

target_link_libraries (NETGeographic Geographic)

# Set the version number on the library
set_target_properties (NETGeographic PROPERTIES VERSION "${LIBVERSIONFULL}" )

# Specify where the library is installed, adding it to the export depends
install (TARGETS NETGeographic
  EXPORT depends
  RUNTIME DESTINATION bin)

if (MSVC)
  # Install pdb files in debug mode
  get_property(PDB TARGET Geographic PROPERTY LOCATION_DEBUG)
  get_filename_component(_P ${PDB} PATH)
  get_filename_component(PDB ${PDB} NAME_WE)
  set(PDB ${_P}/${PDB}.pdb)
  INSTALL(FILES ${PDB} DESTINATION lib${LIB_SUFFIX} CONFIGURATIONS Debug)
endif ()

# Put all the library into a folder in the IDE
set_property (TARGET NETGeographic PROPERTY FOLDER library)
