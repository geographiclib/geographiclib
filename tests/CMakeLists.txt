# Build the test programs...

set (TESTPROGRAMS
  ProjTest TMTest GeodTest ConicTest NaNTester HarmTest EllipticTest intersect)

# Check whether the C++11 random routines are available.
check_cxx_source_compiles (
  "#include <random>
int main() {
  unsigned s = std::random_device()(); // Set seed from random_device
  std::mt19937 g(s);                   // Initialize URNG
  std::uniform_real_distribution<double> dis;
  return int(dis(g));
}\n" CXX11_RANDOM)

if (CXX11_RANDOM)
  set (TESTPROGRAMS ${TESTPROGRAMS} GeodShort NormalTest)
endif ()

if (Boost_FOUND AND NOT GEOGRAPHICLIB_PRECISION EQUAL 4)
  # Skip LevelEllipsoid for quad precision because of compiler errors
  # with boost 1.69 and g++ 9.2.1 (Fedora 30).  Problem reported as
  # https://github.com/boostorg/odeint/issues/40
  set (TESTPROGRAMS ${TESTPROGRAMS} LevelEllipsoid)
  include_directories ("${Boost_INCLUDE_DIRS}")
  if (APPLE)
    # Suppress warnings from Boost library
    # warnings with Mac OS X and boost 1.63
    # no warnings with Linux and boost 1.60
    set (CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -Wno-unused-parameter -Wno-unused-variable -Wno-unused-local-typedef")
  endif ()
endif ()

# Loop over all the tools, specifying the source and library.
add_custom_target (testprograms)
foreach (TESTPROGRAM ${TESTPROGRAMS})

  add_executable (${TESTPROGRAM} EXCLUDE_FROM_ALL ${TESTPROGRAM}.cpp)
  add_dependencies (testprograms ${TESTPROGRAM})
  target_link_libraries (${TESTPROGRAM} ${PROJECT_LIBRARIES}
    ${HIGHPREC_LIBRARIES})

endforeach ()

add_executable (GeodExact EXCLUDE_FROM_ALL GeodExact.cpp
  Geodesic30.cpp GeodesicLine30.cpp
  Geodesic30.hpp GeodesicLine30.hpp)
add_dependencies (testprograms GeodExact)
target_link_libraries (GeodExact ${PROJECT_LIBRARIES} ${HIGHPREC_LIBRARIES})
set (TESTPROGRAMS ${TESTPROGRAMS} GeodExact)

add_executable (reformat EXCLUDE_FROM_ALL reformat.cpp)
if (MSVC OR CMAKE_CONFIGURATION_TYPES)
  # Add _d suffix for your debug versions of the tools
  set_target_properties (GeodExact PROPERTIES
    DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
endif ()
add_dependencies (testprograms reformat)
set (TESTPROGRAMS ${TESTPROGRAMS} reformat)

if (MSVC OR CMAKE_CONFIGURATION_TYPES)
  # Add _d suffix for your debug versions of the tools
  set_target_properties (${TESTPROGRAMS} PROPERTIES
    DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
endif ()

# Put all the tools into a folder in the IDE
set_property (TARGET testprograms ${TESTPROGRAMS} PROPERTY FOLDER tests)

# Don't install test programs
