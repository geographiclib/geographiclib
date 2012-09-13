A library for geographic projections.

Written by Charles Karney <charles@karney.com> and licensed under
the MIT/X11 License.  For more information, see

    http://geographiclib.sourceforge.net/

Files

    00README.txt  -- this file
    AUTHORS -- the authors of the library
    LICENSE.txt -- the MIT/X11 License
    INSTALL -- brief installation instructions
    NEWS -- a history of changes

    include/GeographicLib/ and src/
      Config.h.in, Config.h -- system dependent configuration
      Constants.hpp -- WGS84 constants
      Math.hpp -- math routines
      Accumulator.hpp -- quad precision adder
      Utility.hpp -- I/O and date routines
      PolarStereographic.[ch]pp -- polar stereographic projection
      TransverseMercator.[ch]pp -- transverse Mercator projection
      UTMUPS.[ch]pp -- UTM and UPS
      MGRS.[ch]pp -- MGRS
      TransverseMercatorExact.[ch]pp -- exact TM projection
      EllipticFunction.[ch]pp -- elliptic functions
      GeoCoords.[ch]pp -- hold geographic location
      DMS.[ch]pp -- handle degrees minutes seconds
      Geocentric.[ch]pp -- geocentric coordinates
      LocalCartesian.[ch]pp -- local cartesian coordinates
      Geodesic.[ch]pp -- geodesic calculations
      GeodesicLine.[ch]pp -- calculations on a single geodesic
      PolygonArea.[ch]pp -- polygon area
      AzimuthalEquidistant.[ch]pp -- azimuthal equidistant projection
      Gnomonic.[ch]pp -- gnomonic projection
      CassiniSoldner.[ch]pp -- Cassini-Soldner equidistant projection
      Geoid.[ch]pp -- geoid heights
      Gravity{Model,Circle}.[ch]pp -- gravity models
      Magnetic{Model,Circle}.[ch]pp -- geomagentic models
      {Spherical,Circular}Engine.[ch]pp -- spherical harmonic sums
      SphericalHarmonic{,1,2}.hpp -- frontend for spherical harmonics
      LambertConformalConic.[ch]pp -- Lambert conformal conic projection
      AlbersEqualArea.[ch]pp -- Albers equal area projection
      Gnomonic.[ch]pp -- Ellipsoidal gnomonic projection
      OSGB.[ch]pp -- Ordnance Survey grid system
      Geohash.[ch]pp -- conversions for geohashes
      Ellipsoid.[ch]pp -- ellipsoid properties

    examples/
      example-*.cpp -- simple usage examples for all the classes
      GeoidToGTX.cpp -- a parallelization example

    tools/
      GeoConvert.cpp -- geographic conversion utility
      TransverseMercatorTest.cpp -- TM tester
      Geod.cpp -- geodesic utility
      CartConvert.cpp -- convert to geocentric and local cartesian
      EquidistantTest.cpp -- exercise AzimuthalEquidistant and CassiniSoldner
      GeoidEval.cpp -- evaluate geoid heights
      Gravity.cpp -- evaluate gravity
      MagneticField.cpp -- evaluate magnetic field
      Planimeter.cpp -- computer polygon areas
      geographiclib-get-geoids -- download geoid datasets
      geographiclib-get-magnetic -- download geomagnetic models

    windows/
      GeographicLib-vc9.sln -- MS Studio 2008 solution
      Geographic-vc9.vcproj -- project for library
      GeoConvert-vc9.vcproj -- project for GeoConvert
      TransverseMercatorTest-vc9.vcproj -- project for TransverseMercatorTest
      Geod-vc9.vcproj -- project for Geod
      Planimeter-vc9.vcproj -- project for Planimeter
      CartConvert-vc9.vcproj -- project for CartConvert
      EquidistantTest-vc9.vcproj -- project for EquidistantTest
      GeoidEval-vc9.vcproj -- project for GeoidEval
      Gravity-vc9.vcproj -- project for Gravity
      MagneticField-vc9.vcproj -- project for MagneticField
      also files for MS Studio 2005 (with vc8)
      also files for MS Studio 2010 (with vc10)

    maxima/
      tm.mac -- Maxima code for high precision TM
      ellint.mac -- Maxima code for elliptic functions needed by tm.mac
      tmseries.mac -- Maxima code for series approximations for TM
      geod.mac -- Maxima code for series approximates for Geodesic

    matlab/
      geographiclibinterface.m -- Matlab code to compile Matlab interfaces
      utmupsforward.{cpp,m} -- Matlab code to convert geographic to UTM/UPS
      utmupsreverse.{cpp,m} -- Matlab code to convert UTM/UPS to geographic
      mgrsforward.{cpp,m} -- Matlab code to convert UTM/UPS to MGRS
      mgrsreverse.{cpp,m} -- Matlab code to convert MGRS to UTM/UPS
      geodesicdirect.{cpp,m} -- Matlab code for the direct geodesic problem
      geodesicinverse.{cpp,m} -- Matlab code for the inverse geodesic problem
      geodesicline.{cpp,m} -- Matlab code for geodesic lines
      geoidheight.{cpp,m} -- Matlab code to look up geoid heights
      polygonarea.{cpp,m} -- Matlab code for polygon areas

    doc/
      doxyfile.in -- Doxygen config file
      Geographic.dox -- main page of Doxygen documentation
      geodseries30.html -- geodesic series to 30th order
      tmseries30.html -- transverse Mercator series to 30th order
      html/* -- directory with built documentation
      scripts/*.html -- demonstrations of the Javascript interface
      scripts/GeographicLib/*.js -- Javascript implementation of geodesics

    man/
      *.pod -- plain old documentation
      *.1 -- man pages in nroff format
      *.1.html -- man pages in html format
      *.usage -- documentation for incorporation into executables

    python/GeographicLib/*.py -- Python implementation of geodesics

    Makefile.mk -- Unix/Linux makefiles
    configure -- autoconf configuration script
    CMakeLists.txt -- cmake configuration files
    cmake/
      FindGeographicLib.cmake -- cmake find script
      *.cmake.in -- cmake config templates
