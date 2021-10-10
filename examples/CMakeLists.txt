# Compile a bunch of tiny example programs.  These are built with the
# "exampleprograms" target.  These are mainly for including as examples
# within the doxygen documentation; however, compiling them catches some
# obvious blunders.
if (GEOGRAPHICLIB_PRECISION EQUAL 2)
  # These examples all assume real = double
  set (EXAMPLE_SOURCES
    example-Accumulator.cpp
    example-AlbersEqualArea.cpp
    example-AzimuthalEquidistant.cpp
    example-CassiniSoldner.cpp
    example-CircularEngine.cpp
    example-Constants.cpp
    example-DMS.cpp
    example-Ellipsoid.cpp
    example-EllipticFunction.cpp
    example-GARS.cpp
    example-GeoCoords.cpp
    example-Geocentric.cpp
    example-Geodesic.cpp
    example-Geodesic-small.cpp
    example-GeodesicExact.cpp
    example-GeodesicLine.cpp
    example-GeodesicLineExact.cpp
    example-GeographicErr.cpp
    example-Geohash.cpp
    example-Geoid.cpp
    example-Georef.cpp
    example-Gnomonic.cpp
    example-GravityCircle.cpp
    example-GravityModel.cpp
    example-LambertConformalConic.cpp
    example-LocalCartesian.cpp
    example-MGRS.cpp
    example-MagneticCircle.cpp
    example-MagneticModel.cpp
    example-Math.cpp
    example-NearestNeighbor.cpp
    example-NormalGravity.cpp
    example-OSGB.cpp
    example-PolarStereographic.cpp
    example-PolygonArea.cpp
    example-Rhumb.cpp
    example-RhumbLine.cpp
    example-SphericalEngine.cpp
    example-SphericalHarmonic.cpp
    example-SphericalHarmonic1.cpp
    example-SphericalHarmonic2.cpp
    example-TransverseMercator.cpp
    example-TransverseMercatorExact.cpp
    example-UTMUPS.cpp
    example-Utility.cpp
    )
  if (USE_BOOST_FOR_EXAMPLES AND Boost_FOUND)
    add_definitions (-DGEOGRAPHICLIB_HAVE_BOOST_SERIALIZATION=1)
    include_directories ("${Boost_INCLUDE_DIRS}")
  endif ()
else ()
  set (EXAMPLE_SOURCES)
endif ()
set (EXAMPLE_SOURCES ${EXAMPLE_SOURCES}
  GeoidToGTX.cpp make-egmcof.cpp JacobiConformal.cpp)

set (EXAMPLES)
add_definitions (${PROJECT_DEFINITIONS})

foreach (EXAMPLE_SOURCE ${EXAMPLE_SOURCES})
  get_filename_component (EXAMPLE ${EXAMPLE_SOURCE} NAME_WE)
  set (EXAMPLES ${EXAMPLES} ${EXAMPLE})
  if (EXISTS ${EXAMPLE}.hpp)
    set (EXAMPLE_SOURCE ${EXAMPLE_SOURCE} ${EXAMPLE}.hpp)
  endif ()
  add_executable (${EXAMPLE} EXCLUDE_FROM_ALL ${EXAMPLE_SOURCE})
  target_link_libraries (${EXAMPLE} ${PROJECT_LIBRARIES} ${HIGHPREC_LIBRARIES})
endforeach ()

if (Boost_FOUND AND GEOGRAPHICLIB_PRECISION EQUAL 2)
  target_link_libraries (example-NearestNeighbor ${Boost_LIBRARIES})
endif ()

find_package (OpenMP QUIET)

if (OPENMP_FOUND OR OpenMP_FOUND)
  set_target_properties (GeoidToGTX PROPERTIES
    COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
  if (NOT WIN32)
    set_target_properties (GeoidToGTX PROPERTIES
      LINK_FLAGS ${OpenMP_CXX_FLAGS})
  endif ()
endif ()

if (MSVC OR CMAKE_CONFIGURATION_TYPES)
  # Add _d suffix for your debug versions of the tools
  set_target_properties (${EXAMPLES} PROPERTIES
    DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
endif ()

add_custom_target (exampleprograms DEPENDS ${EXAMPLES})

# Put all the examples into a folder in the IDE
set_property (TARGET exampleprograms ${EXAMPLES} PROPERTY FOLDER examples)
