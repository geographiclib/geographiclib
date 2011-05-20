# $Id: 00README.txt 6867 2010-09-11 13:04:26Z karney $

A library for geographic projections.

Written by Charles Karney <charles@karney.com> and licensed under
the LGPL.  For more information, see

    http://geographiclib.sourceforge.net/

Files

    00README.txt  -- this file
    AUTHORS -- the authors of the library
    COPYING -- the LGPL license, v. 3
    INSTALL -- brief installation instructions
    NEWS -- a history of changes

    include/GeographicLib/ and src/
      Constants.hpp -- WGS84 constants
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
      AzimuthalEquidistant.[ch]pp -- azimuthal equidistant projection
      Gnomonic.[ch]pp -- gnomonic projection
      CassiniSoldner.[ch]pp -- Cassini-Soldner equidistant projection
      Geoid.[ch]pp -- geoid heights
      LambertConformalConic.[ch]pp -- Lambert conformal conic projection
      Gnomonic.[ch]pp -- Ellipsoidal gnomonic projection

    tools/
      GeoConvert.cpp -- geographic conversion utility
      TransverseMercatorTest.cpp -- TM tester
      Geod.cpp -- geodesic utility
      CartConvert.cpp -- convert to geocentric and local cartesian
      EquidistantTest.cpp -- exercise AzimuthalEquidistant and CassiniSoldner
      GeoidEval.cpp -- evaluate geoid heights
      Planimeter.cpp -- computer polygon areas

    windows/
      GeographicLib-vc9.sln -- MS Studio 2008 solution
      GeographicLib-vc9.vcproj -- project for library
      GeoConvert-vc9.vcproj -- project for GeoConvert
      TransverseMercatorTest-vc9.vcproj -- project for TransverseMercatorTest
      Geod-vc9.vcproj -- project for Geod
      Planimeter-vc9.vcproj -- project for Planimeter
      CartConvert-vc9.vcproj -- project for CartConvert
      EquidistantTest-vc9.vcproj -- project for EquidistantTest
      GeoidEval-vc9.vcproj -- project for GeoidEval
      also files for MS Studio 2005 (with vc8)

    maxima/
      tm.mac -- Maxima code for high precision TM
      ellint.mac -- Maxima code for elliptic functions needed by tm.mac
      tmseries.mac -- Maxima code for series approximations for TM
      geod.mac -- Maxima code for series approximates for Geodesic

    doc/
      Doxyfile -- Doxygen config file
      Geographic.doc -- main page of Doxygen documentation
      gauss-laborde-graticule-a.{png,pdf} -- Fig. 1
      gauss-krueger-graticule-a.{png,pdf} -- Fig. 2
      thompson-tm-graticule-a.{png,pdf} -- Fig. 3
      gauss-krueger-graticule.{png,pdf} -- Fig. 4
      gauss-krueger-convergence-scale.{png,pdf} -- Fig. 5
      thompson-tm-graticule.{png,pdf} -- Fig. 6
      html/* -- directory with built documentation

    Makefile -- Unix/Linux makefile (invokes make in the other directories)
    configure -- autoconf configuration script
