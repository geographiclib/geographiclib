# $Id$

A library for geographic projections.

Written by Charles Karney <charles@karney.com> and licensed under
the LGPL.  For more information, see

    http://charles.karney.info/geographic/

Files

    00README.txt  -- this file
    Doxyfile -- Doxygen config file
    Geographic.doc -- main page of Doxygen documentation

    PolarStereographic.[ch]pp -- polar stereographic projection
    TransverseMercator.[ch]pp -- transverse Mercator projection
    UTMUPS.[ch]pp -- UTM and UPS
    MGRS.[ch]pp -- MGRS
    Constants.[ch]pp -- WGS84 constants
    TransverseMercatorExact.[ch]pp -- exact TM projection
    EllipticFunction.[ch]pp -- elliptic functions
    GeoCoords.[ch]pp -- hold geographic location
    DMS.[ch]pp -- handle degrees minutes seconds
    Geocentric.[ch]pp -- geocentric coordinates
    LocalCartesian.[ch]pp -- local cartesian coordinates
    Geodesic.[ch]pp -- geodesic calculatiosn

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

    gauss-laborde-graticule-a.{png,pdf} -- Fig. 1
    gauss-krueger-graticule-a.{png,pdf} -- Fig. 2
    thompson-tm-graticule-a.{png,pdf} -- Fig. 3
    gauss-krueger-graticule.{png,pdf} -- Fig. 4
    gauss-krueger-convergence-scale.{png,pdf} -- Fig. 5
    thompson-tm-graticule.{png,pdf} -- Fig. 6

This is the 2009-03 version of the library.

Changes between 2009-03 and 2009-02 versions:

  * Rename Geographic::ECEF to Geographic::Geocentric, ECEFConvert to
    CartConvert

  * Add Geographic::Geodesic and the Geod utility.

  * Declare when no exceptions are thrown by functions.

  * Minor changes to GeographicLib::DMS class.

  * Use invf = 0 to mean a sphere in constructors to some classes.

  * The makefile creates a library and includes an install target.

Changes between 2009-02 and 2009-01 versions:

  * Fix documentation of constructors (flattening -> inverse
    flattening).

  * Use std versions of math functions.

  * Add ECEF and LocalCartesian classes and ECEFConvert utility.

  * Gather the documentation on the utility programs onto one page.
