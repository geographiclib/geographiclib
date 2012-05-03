# $Id$

MODULES = DMS EllipticFunction GeoCoords MGRS PolarStereographic \
	TransverseMercator TransverseMercatorExact UTMUPS Geocentric \
	LocalCartesian Geodesic GeodesicLine PolygonArea \
	AzimuthalEquidistant CassiniSoldner \
	Geoid Gnomonic OSGB AlbersEqualArea Geohash
PROGRAMS = GeoConvert TransverseMercatorProj CartConvert Geod GeodesicProj \
	GeoidEval Gravity MagneticField Planimeter ConicProj

HEADERS = $(patsubst %,../include/GeographicLib/%.hpp,Constants $(MODULES))
SOURCES = $(patsubst %,../src/%.cpp,$(MODULES)) \
	$(patsubst %,../tools/%.cpp,$(PROGRAMS))

EXTRAFILES = tmseries30.html geodseries30.html
HTMLMANPAGES = 	$(patsubst %,../man/%.1.html,$(PROGRAMS))
SCRIPTDRIVERS = $(wildcard scripts/*.html)
JSSCRIPTS = $(wildcard scripts/GeographicLib/*.js)

MAXIMA = tm ellint tmseries geod
MAXIMASOURCES = $(patsubst %,../maxima/%.mac,$(MAXIMA))

VERSION:=$(shell grep '\bVERSION=' ../configure | cut -f2 -d\' | head -1)

doc: html/index.html

html/index.html: doxyfile.in Geographic.doc \
	$(HEADERS) $(ALLSOURCES) $(MAXIMASOURCES) $(EXTRAFILES) \
	$(HTMLMANPAGES)
	if test -d html; then rm -rf html/*; else mkdir html; fi
	cp -p $(MAXIMASOURCES) $(EXTRAFILES) $(HTMLMANPAGES) \
	../LICENSE.txt html/
	sed -e "s%@PROJECT_SOURCE_DIR@%..%g" \
	-e "s%@GeographicLib_VERSION@%$(VERSION)%g" \
	doxyfile.in | doxygen -

PREFIX = /usr/local
DEST = $(PREFIX)/share/doc/GeographicLib
DOCDEST = $(DEST)/html
SCRIPTDEST = $(DEST)/scripts
INSTALL = install -b

install: html/index.html
	test -d $(DOCDEST) || mkdir -p $(DOCDEST)
	$(INSTALL) -m 644 html/* $(DOCDEST)/
	test -d $(SCRIPTDEST)/GeographicLib || \
	mkdir -p $(SCRIPTDEST)/GeographicLib
	$(INSTALL) -m 644 $(SCRIPTDRIVERS) $(SCRIPTDEST)/
	$(INSTALL) -m 644 $(JSSCRIPTS) $(SCRIPTDEST)/GeographicLib/

maintainer-clean:
	rm -rf html

.PHONY: doc install clean
