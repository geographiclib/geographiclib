# $Id: Makefile.mk 6906 2010-12-02 22:10:56Z karney $

MODULES = DMS EllipticFunction GeoCoords MGRS PolarStereographic \
	TransverseMercator TransverseMercatorExact UTMUPS Geocentric \
	LocalCartesian Geodesic GeodesicLine \
	AzimuthalEquidistant CassiniSoldner \
	Geoid LambertConformalConic Gnomonic OSGB AlbersEqualArea

PREFIX = /usr/local
LIBNAME = GeographicLib
HEADERS = $(LIBNAME)/Constants.hpp $(patsubst %,$(LIBNAME)/%.hpp,$(MODULES))
DEST = $(PREFIX)/include/$(LIBNAME)

INSTALL=install -b

all:
	@:

install:
	test -d $(DEST) || mkdir -p $(DEST)
	$(INSTALL) -m 644 $(HEADERS) $(DEST)/
list:
	@echo $(HEADERS)

clean:
	@:

.PHONY: all install list clean
