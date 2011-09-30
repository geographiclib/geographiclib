# $Id: da35a9a983cd00affd326d1364381995512affde $

MODULES = DMS EllipticFunction GeoCoords MGRS PolarStereographic \
	TransverseMercator TransverseMercatorExact UTMUPS Geocentric \
	LocalCartesian Geodesic GeodesicLine PolygonArea \
	AzimuthalEquidistant CassiniSoldner \
	Geoid LambertConformalConic Gnomonic OSGB AlbersEqualArea
EXTRAHEADERS = Constants Math Accumulator

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
