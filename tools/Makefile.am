#
# Makefile.am
#
# Copyright (C) 2009, Francesco P. Lovergine <frankie@debian.org>
# $Id$

AM_CPPFLAGS = -Wall -funroll-loops -finline-functions -fomit-frame-pointer
INCLUDES = -I$(top_srcdir)/include -I$(top_srcdir)/man

LDADD = $(top_builddir)/src/libGeographic.la
DEPS = $(top_builddir)/src/libGeographic.la

bin_PROGRAMS = GeoConvert \
	TransverseMercatorProj \
	CartConvert \
	Geod \
	GeodesicProj \
	GeoidEval \
	Planimeter \
	ConicProj

GeoConvert.o: ../man/GeoConvert.usage
TransverseMercatorProj.o: ../man/TransverseMercatorProj.usage
CartConvert.o: ../man/CartConvert.usage
Geod.o: ../man/Geod.usage
GeodesicProj.o: ../man/GeodesicProj.usage
GeoidEval.o: ../man/GeoidEval.usage
Planimeter.o: ../man/Planimeter.usage
ConicProj.o: ../man/ConicProj.usage

GeoConvert_SOURCES = GeoConvert.cpp \
	../include/GeographicLib/Constants.hpp \
	../include/GeographicLib/Math.hpp \
	../include/GeographicLib/DMS.hpp \
	../include/GeographicLib/GeoCoords.hpp \
	../include/GeographicLib/UTMUPS.hpp

TransverseMercatorProj_SOURCES = TransverseMercatorProj.cpp \
	../include/GeographicLib/Constants.hpp \
	../include/GeographicLib/Math.hpp \
	../include/GeographicLib/DMS.hpp \
	../include/GeographicLib/EllipticFunction.hpp \
	../include/GeographicLib/TransverseMercator.hpp \
	../include/GeographicLib/TransverseMercatorExact.hpp

CartConvert_SOURCES = CartConvert.cpp \
	../include/GeographicLib/Constants.hpp \
	../include/GeographicLib/Math.hpp \
	../include/GeographicLib/DMS.hpp \
	../include/GeographicLib/Geocentric.hpp \
	../include/GeographicLib/LocalCartesian.hpp

Geod_SOURCES = Geod.cpp \
	../include/GeographicLib/Constants.hpp \
	../include/GeographicLib/Math.hpp \
	../include/GeographicLib/DMS.hpp \
	../include/GeographicLib/Geodesic.hpp \
	../include/GeographicLib/GeodesicLine.hpp

GeodesicProj_SOURCES = GeodesicProj.cpp \
	../include/GeographicLib/AzimuthalEquidistant.hpp \
	../include/GeographicLib/CassiniSoldner.hpp \
	../include/GeographicLib/Gnomonic.hpp \
	../include/GeographicLib/Constants.hpp \
	../include/GeographicLib/Math.hpp \
	../include/GeographicLib/DMS.hpp \
	../include/GeographicLib/Geodesic.hpp \
	../include/GeographicLib/GeodesicLine.hpp

GeoidEval_SOURCES = GeoidEval.cpp \
	../include/GeographicLib/Constants.hpp \
	../include/GeographicLib/Math.hpp \
	../include/GeographicLib/DMS.hpp \
	../include/GeographicLib/GeoCoords.hpp \
	../include/GeographicLib/Geoid.hpp

Planimeter_SOURCES = Planimeter.cpp \
	../include/GeographicLib/Constants.hpp \
	../include/GeographicLib/Math.hpp \
	../include/GeographicLib/DMS.hpp \
	../include/GeographicLib/GeoCoords.hpp \
	../include/GeographicLib/PolygonArea.hpp \
	../include/GeographicLib/Accumulator.hpp \
	../include/GeographicLib/Geodesic.hpp

ConicProj_SOURCES = ConicProj.cpp \
	../include/GeographicLib/Constants.hpp \
	../include/GeographicLib/Math.hpp \
	../include/GeographicLib/DMS.hpp \
	../include/GeographicLib/AlbersEqualArea.hpp \
	../include/GeographicLib/LambertConformalConic.hpp

sbin_SCRIPTS = geographiclib-get-geoids

geoids_dir = $(datadir)/GeographicLib/geoids
geographiclib-get-geoids: geographiclib-get-geoids.sh
	sed -e "s%@GEOID_DEFAULT_PATH@%$(geoids_dir)%" $< > $@
	chmod +x $@

CLEANFILES = geographiclib-get-geoids
EXTRA_DIST = Makefile.mk CMakeLists.txt geographiclib-get-geoids.sh
