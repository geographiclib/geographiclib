# $Id: 73c53dea178abd7300904ae145208ead071960de $

MODULES = DMS EllipticFunction GeoCoords MGRS PolarStereographic \
	TransverseMercator TransverseMercatorExact UTMUPS Geocentric \
	LocalCartesian Geodesic GeodesicLine \
	AzimuthalEquidistant CassiniSoldner \
	Geoid LambertConformalConic Gnomonic OSGB AlbersEqualArea

PREFIX = /usr/local
LIBNAME = GeographicLib
HEADERS = $(LIBNAME)/Constants.hpp $(LIBNAME)/Config.h \
	$(patsubst %,$(LIBNAME)/%.hpp,$(MODULES))
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
