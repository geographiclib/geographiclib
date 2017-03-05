# Build the test programs...

set (TESTPROGRAMS
  ProjTest TMTest GeodTest ConicTest NaNTester HarmTest EllipticTest)

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

if (Boost_FOUND)
  set (TESTPROGRAMS ${TESTPROGRAMS} LevelEllipsoid)
  include_directories ("${Boost_INCLUDE_DIRS}")
  # Suppress warnings from Boost library
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter -Wno-unused-variable -Wno-unused-local-typedef")
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

# Put all the tools into a folder in the IDE
set_property (TARGET testprograms ${TESTPROGRAMS} PROPERTY FOLDER tests)

# Don't install test programs
