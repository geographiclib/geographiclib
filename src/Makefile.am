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
		DMS.cpp \
		EllipticFunction.cpp \
		Geocentric.cpp \
		GeoCoords.cpp \
		Geodesic.cpp \
		GeodesicLine.cpp \
		Geoid.cpp \
		Gnomonic.cpp \
		LambertConformalConic.cpp \
		LocalCartesian.cpp \
		MGRS.cpp \
		OSGB.cpp \
		PolarStereographic.cpp \
		PolygonArea.cpp \
		TransverseMercator.cpp \
		TransverseMercatorExact.cpp \
		UTMUPS.cpp \
		../include/GeographicLib/Accumulator.hpp \
		../include/GeographicLib/AlbersEqualArea.hpp \
		../include/GeographicLib/AzimuthalEquidistant.hpp \
		../include/GeographicLib/CassiniSoldner.hpp \
		../include/GeographicLib/Constants.hpp \
		../include/GeographicLib/DMS.hpp \
		../include/GeographicLib/EllipticFunction.hpp \
		../include/GeographicLib/Geocentric.hpp \
		../include/GeographicLib/GeoCoords.hpp \
		../include/GeographicLib/Geodesic.hpp \
		../include/GeographicLib/GeodesicLine.hpp \
		../include/GeographicLib/Geoid.hpp \
		../include/GeographicLib/Gnomonic.hpp \
		../include/GeographicLib/LambertConformalConic.hpp \
		../include/GeographicLib/LocalCartesian.hpp \
		../include/GeographicLib/Math.hpp \
		../include/GeographicLib/MGRS.hpp \
		../include/GeographicLib/OSGB.hpp \
		../include/GeographicLib/PolarStereographic.hpp \
		../include/GeographicLib/PolygonArea.hpp \
		../include/GeographicLib/TransverseMercatorExact.hpp \
		../include/GeographicLib/TransverseMercator.hpp \
		../include/GeographicLib/UTMUPS.hpp \
		../include/GeographicLib/Config.h

../include/GeographicLib/Config.h: ../include/GeographicLib/Config-ac.h
	egrep '\bVERSION\b|\bHAVE_LONG_DOUBLE\b' $< | \
	sed -e s/VERSION/GEOGRAPHICLIB_VERSION_STRING/ > $@

$(libGeographic_la_OBJECTS): ../include/GeographicLib/Config.h

geoids_dir=$(datadir)/GeographicLib/geoids

DEFS=-DGEOID_DEFAULT_PATH=\"$(geoids_dir)\"  @DEFS@

EXTRA_DIST = Makefile.mk CMakeLists.txt
