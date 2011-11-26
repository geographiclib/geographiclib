#
# Makefile.am
#
# Copyright (C) 2009, Francesco P. Lovergine <frankie@debian.org>
# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/include -Wall -Wextra

lib_LTLIBRARIES = libGeographic.la

libGeographic_la_LDFLAGS = -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)
libGeographic_la_SOURCES = AlbersEqualArea.cpp \
		AzimuthalEquidistant.cpp \
		CassiniSoldner.cpp \
		CircularEngine.cpp \
		DMS.cpp \
		EllipticFunction.cpp \
		GeoCoords.cpp \
		Geocentric.cpp \
		Geodesic.cpp \
		GeodesicLine.cpp \
		Geoid.cpp \
		Gnomonic.cpp \
		LambertConformalConic.cpp \
		LocalCartesian.cpp \
		MGRS.cpp \
		MagneticCircle.cpp \
		MagneticModel.cpp \
		OSGB.cpp \
		PolarStereographic.cpp \
		PolygonArea.cpp \
		SphericalEngine.cpp \
		TransverseMercator.cpp \
		TransverseMercatorExact.cpp \
		UTMUPS.cpp \
		../include/GeographicLib/Accumulator.hpp \
		../include/GeographicLib/AlbersEqualArea.hpp \
		../include/GeographicLib/AzimuthalEquidistant.hpp \
		../include/GeographicLib/CassiniSoldner.hpp \
		../include/GeographicLib/CircularEngine.hpp \
		../include/GeographicLib/Constants.hpp \
		../include/GeographicLib/DMS.hpp \
		../include/GeographicLib/EllipticFunction.hpp \
		../include/GeographicLib/GeoCoords.hpp \
		../include/GeographicLib/Geocentric.hpp \
		../include/GeographicLib/Geodesic.hpp \
		../include/GeographicLib/GeodesicLine.hpp \
		../include/GeographicLib/Geoid.hpp \
		../include/GeographicLib/Gnomonic.hpp \
		../include/GeographicLib/LambertConformalConic.hpp \
		../include/GeographicLib/LocalCartesian.hpp \
		../include/GeographicLib/MGRS.hpp \
		../include/GeographicLib/MagneticCircle.hpp \
		../include/GeographicLib/MagneticModel.hpp \
		../include/GeographicLib/Math.hpp \
		../include/GeographicLib/OSGB.hpp \
		../include/GeographicLib/PolarStereographic.hpp \
		../include/GeographicLib/PolygonArea.hpp \
		../include/GeographicLib/SphericalEngine.hpp \
		../include/GeographicLib/SphericalHarmonic.hpp \
		../include/GeographicLib/SphericalHarmonic1.hpp \
		../include/GeographicLib/SphericalHarmonic2.hpp \
		../include/GeographicLib/TransverseMercator.hpp \
		../include/GeographicLib/TransverseMercatorExact.hpp \
		../include/GeographicLib/UTMUPS.hpp \
		../include/GeographicLib/Utility.hpp \
		../include/GeographicLib/Config.h

../include/GeographicLib/Config.h: ../include/GeographicLib/Config-ac.h
	( egrep '\bVERSION\b|\bHAVE_LONG_DOUBLE\b' $< | \
	  sed -e s/VERSION/GEOGRAPHICLIB_VERSION_STRING/; \
	  grep WORDS_BIGENDIAN $< | tail -1 ) > $@

$(libGeographic_la_OBJECTS): ../include/GeographicLib/Config.h

geographiclib_data=$(datadir)/GeographicLib

DEFS=-DGEOGRAPHICLIB_DATA=\"$(geographiclib_data)\" @DEFS@

EXTRA_DIST = Makefile.mk CMakeLists.txt Geographic.pro
