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
    ECEF.[ch]pp -- ECEF coordinates
    LocalCartesian.[ch]pp -- local cartesian coordinates

    GeoConvert.cpp -- geographic conversion utility
    TransverseMercatorTest.cpp -- TM tester
    ECEFConvert.cpp -- convert to ECEF and local cartesian

    Makefile -- Unix/Linux makefile

    GeographicLib.sln -- MS Studio 2005 solution
    GeographicLib.vcproj -- project for library
    GeoConvert.vcproj -- project for GeoConvert
    TransverseMercatorTest.vcproj -- project for TransverseMercatorTest
    ECEFConvert.vcproj -- project for ECEFConvert

    tm.mac -- Maxima code for high precision TM
    ellint.mac -- Maxima code for elliptic functions needed by tm.mac
    tmseries.mac -- Maxima code for series approximations for TM

    gauss-laborde-graticule-a.{png,pdf} -- Fig. 1
    gauss-krueger-graticule-a.{png,pdf} -- Fig. 2
    thompson-tm-graticule-a.{png,pdf} -- Fig. 3
    gauss-krueger-graticule.{png,pdf} -- Fig. 4
    gauss-krueger-convergence-scale.{png,pdf} -- Fig. 5
    thompson-tm-graticule.{png,pdf} -- Fig. 6

This is the 2009-02 version of the library.

Changes between 2009-02 and 2009-01 versions:

  * Fix documentation of constructors (flattening -> inverse
    flattening).

  * Use std versions of math functions.

  * Add ECEF and LocalCartesian classes and ECEFConvert utility.

  * Gather the documentation on the utility programs onto one page.
