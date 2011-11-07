# $Id$

MODULES = AlbersEqualArea AzimuthalEquidistant CassiniSoldner \
	CircularEngine DMS EllipticFunction GeoCoords Geocentric \
	Geodesic GeodesicLine Geoid Gnomonic LambertConformalConic \
	LocalCartesian MGRS MagneticCircle MagneticModel OSGB \
	PolarStereographic PolygonArea SphericalEngine \
	TransverseMercator TransverseMercatorExact \
	UTMUPS
EXTRAHEADERS = Accumulator Constants Math SphericalHarmonic \
	SphericalHarmonic1 SphericalHarmonic2 Utility

PREFIX = /usr/local
LIBNAME = GeographicLib
HEADERS = $(LIBNAME)/Config.h \
	$(patsubst %,$(LIBNAME)/%.hpp,$(EXTRAHEADERS) $(MODULES))
DEST = $(PREFIX)/include/$(LIBNAME)

INSTALL = install -b

all:
	@:

install:
	test -d $(DEST) || mkdir -p $(DEST)
	$(INSTALL) -m 644 $(HEADERS) $(DEST)/
clean:
	@:

.PHONY: all install clean
