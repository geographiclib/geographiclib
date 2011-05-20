#
# Makefile.am
#
# Copyright (C) 2009, Francesco P. Lovergine <frankie@debian.org>
# $Id: Makefile.am 6861 2010-09-06 21:29:08Z karney $

AM_CPPFLAGS = -Wall -funroll-loops -finline-functions -fomit-frame-pointer
INCLUDES = -I$(top_srcdir)/include 

LDADD = $(top_builddir)/src/libGeographic.la
DEPS = $(top_builddir)/src/libGeographic.la

bin_PROGRAMS = GeoConvert TransverseMercatorTest CartConvert Geod EquidistantTest GeoidEval Planimeter

GeoConvert_SOURCES = GeoConvert.cpp \
	../include/GeographicLib/Constants.hpp \
	../include/GeographicLib/DMS.hpp \
	../include/GeographicLib/GeoCoords.hpp \
	../include/GeographicLib/UTMUPS.hpp

TransverseMercatorTest_SOURCES = TransverseMercatorTest.cpp \
	../include/GeographicLib/Constants.hpp \
	../include/GeographicLib/DMS.hpp \
	../include/GeographicLib/EllipticFunction.hpp \
	../include/GeographicLib/TransverseMercator.hpp \
	../include/GeographicLib/TransverseMercatorExact.hpp

CartConvert_SOURCES = CartConvert.cpp \
	../include/GeographicLib/Constants.hpp \
	../include/GeographicLib/DMS.hpp \
	../include/GeographicLib/Geocentric.hpp \
	../include/GeographicLib/LocalCartesian.hpp

Geod_SOURCES = Geod.cpp \
	../include/GeographicLib/Constants.hpp \
	../include/GeographicLib/DMS.hpp \
	../include/GeographicLib/Geodesic.hpp \
	../include/GeographicLib/GeodesicLine.hpp

EquidistantTest_SOURCES = EquidistantTest.cpp \
	../include/GeographicLib/AzimuthalEquidistant.hpp \
	../include/GeographicLib/CassiniSoldner.hpp \
	../include/GeographicLib/Gnomonic.hpp \
	../include/GeographicLib/Constants.hpp \
	../include/GeographicLib/DMS.hpp \
	../include/GeographicLib/Geodesic.hpp \
	../include/GeographicLib/GeodesicLine.hpp

GeoidEval_SOURCES = GeoidEval.cpp \
	../include/GeographicLib/Constants.hpp \
	../include/GeographicLib/DMS.hpp \
	../include/GeographicLib/GeoCoords.hpp \
	../include/GeographicLib/Geoid.hpp

Planimeter_SOURCES = Planimeter.cpp \
	../include/GeographicLib/Constants.hpp \
	../include/GeographicLib/DMS.hpp \
	../include/GeographicLib/GeoCoords.hpp \
	../include/GeographicLib/Geodesic.hpp \
	../include/GeographicLib/GeodesicLine.hpp

EXTRA_DIST = Makefile.mk
