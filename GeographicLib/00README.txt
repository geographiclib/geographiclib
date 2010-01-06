# $Id$

A library for geographic projections.

Written by Charles Karney <charles@karney.com> and licensed under
the LGPL.  For more information, see

    http://geographiclib.sourceforge.net/

Files

    00README.txt  -- this file
    AUTHORS -- the authors of the library
    NEWS -- a history of changes
    COPYING -- the LGPL license, v. 3

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
      Geodesic.[ch]pp -- geodesic calculatiosn
      AzimuthalEquidistant.[ch]pp -- azimuthal equidistant projection
      CassiniSoldner.[ch]pp -- Cassini-Soldner equidistant projection
      Geoid.[ch]pp -- geoid heights
      LambertConformalConic.[ch]pp -- Lambert conformal conic projection

    tools/
      GeoConvert.cpp -- geographic conversion utility
      TransverseMercatorTest.cpp -- TM tester
      Geod -- geodesic utility
      CartConvert.cpp -- convert to geocentric and local cartesian
      EquidistantTest.cpp -- exercise AzimuthalEquidistant and CassiniSoldner
      GeoidEval.cpp -- evaluate geoid heights

    Makefile -- Unix/Linux makefile (invoke make in the other directories)

    windows/
      GeographicLib.sln -- MS Studio 2005 solution
      GeographicLib.vcproj -- project for library
      GeoConvert.vcproj -- project for GeoConvert
      TransverseMercatorTest.vcproj -- project for TransverseMercatorTest
      Geod.vcproj -- project for Geod
      CartConvert.vcproj -- project for CartConvert
      EquidistantTest.vcproj -- project for EquidistantTest
      GeoidEval.vcproj -- project for GeoidEval

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
