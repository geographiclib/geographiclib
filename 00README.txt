# $Id$

A library for geographic projections.

Written by Charles Karney <charles@karney.com> and licensed under
the LGPL.  For more information, see

    http://charles.karney.info/geographic/

Files

    00README.txt  -- this file
    Doxyfile -- Doxygen config file
    Geographic.doc -- main page of Doxygen documentation

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
    AzimuthalEquidistant.[ch]pp -- azimuthal equidistant coordinates

    GeoConvert.cpp -- geographic conversion utility
    TransverseMercatorTest.cpp -- TM tester
    Geod -- geodesic utility
    CartConvert.cpp -- convert to geocentric and local cartesian

    Makefile -- Unix/Linux makefile

    GeographicLib.sln -- MS Studio 2005 solution
    GeographicLib.vcproj -- project for library
    GeoConvert.vcproj -- project for GeoConvert
    TransverseMercatorTest.vcproj -- project for TransverseMercatorTest
    Geod.vcproj -- project for Geod
    CartConvert.vcproj -- project for CartConvert

    tm.mac -- Maxima code for high precision TM
    ellint.mac -- Maxima code for elliptic functions needed by tm.mac
    tmseries.mac -- Maxima code for series approximations for TM
    geod.mac -- Maxima code for series approximates for Geodesic

    gauss-laborde-graticule-a.{png,pdf} -- Fig. 1
    gauss-krueger-graticule-a.{png,pdf} -- Fig. 2
    thompson-tm-graticule-a.{png,pdf} -- Fig. 3
    gauss-krueger-graticule.{png,pdf} -- Fig. 4
    gauss-krueger-convergence-scale.{png,pdf} -- Fig. 5
    thompson-tm-graticule.{png,pdf} -- Fig. 6

This is the 2009-06 version of the library.

Changes between 2009-06 and 2009-05 versions:

  * Add routines to decode and encode zone+hemisphere to GeographicLib::UTMUPS.

  * Clean up code in GeographicLib::Geodesic.

Changes between 2009-05 and 2009-04 versions:

  * Improvements to GeographicLib::Geodesic:
    + more economical series expansions,
    + return reduced length (as does the Geod utility),
    + improved calculation of starting point for inverse method,
    + use reduced length to give derivative for Newton's method.

  * Add AzimuthalEquidistant class.

  + Make GeographicLib::Geocentric, GeographicLib::TransverseMercator,
    and GeographicLib::PolarStereographic classes work with prolate
    ellipsoids.

  * CartConvert checks its inputs more carefully.

  * Remove reference to defunct Constants.cpp from GeographicLib.vcproj.

Changes between 2009-04 and 2009-03 versions:

  * Use compile-time constants to select the order of series in
    GeographicLib::TransverseMercator.

  * 2x unroll of Clenshaw summation to avoid data shuffling.

  * Simplification of GeographicLib::EllipticFunction::E.

  * Use STATIC_ASSERT for compile-time checking of constants.

  * Improvements to GeographicLib::Geodesic:
    + compile-time option to change order of series used,
    + post maxima code for generating the series,
    + tune the order of series for double,
    + improvements in the selection of starting points for Newton's
      method,
    + accept and return spherical arc lengths,
    + works with both oblate and prolate spheroids,
    + add -a, -e, -b options to the Geod utility.

Changes between 2009-03 and 2009-02 versions:

  * Add GeographicLib::Geodesic and the Geod utility.

  * Declare when no exceptions are thrown by functions.

  * Minor changes to GeographicLib::DMS class.

  * Use invf = 0 to mean a sphere in constructors to some classes.

  * The makefile creates a library and includes an install target.

  * Rename GeographicLib::ECEF to GeographicLib::Geocentric, ECEFConvert
    to CartConvert.

  * Use inline functions to define constant doubles in Constants.hpp.

Changes between 2009-02 and 2009-01 versions:

  * Fix documentation of constructors (flattening -> inverse
    flattening).

  * Use std versions of math functions.

  * Add ECEF and LocalCartesian classes and ECEFConvert utility.

  * Gather the documentation on the utility programs onto one page.
