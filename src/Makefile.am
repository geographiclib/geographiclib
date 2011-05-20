#
# Makefile.am
#
# Copyright (C) 2009, Francesco P. Lovergine <frankie@debian.org>
# $Id: Makefile.am 6859 2010-09-06 14:45:33Z karney $

AM_CPPFLAGS = -I$(top_srcdir)/include -Wall -funroll-loops -finline-functions -fomit-frame-pointer

lib_LTLIBRARIES = libGeographic.la

libGeographic_la_LDFLAGS = -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)	
libGeographic_la_SOURCES = AzimuthalEquidistant.cpp \
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
		PolarStereographic.cpp \
		TransverseMercator.cpp \
		TransverseMercatorExact.cpp \
		UTMUPS.cpp \
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
		../include/GeographicLib/MGRS.hpp \
		../include/GeographicLib/PolarStereographic.hpp \
		../include/GeographicLib/TransverseMercatorExact.hpp \
		../include/GeographicLib/TransverseMercator.hpp \
		../include/GeographicLib/UTMUPS.hpp

geoids_dir=$(datadir)/GeographicLib/geoids

DEFS=-DGEOID_DEFAULT_PATH=\"$(geoids_dir)\"  @DEFS@

EXTRA_DIST = Makefile.mk
