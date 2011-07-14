# $Id$

MODULES = DMS EllipticFunction GeoCoords MGRS PolarStereographic \
	TransverseMercator TransverseMercatorExact UTMUPS Geocentric \
	LocalCartesian Geodesic GeodesicLine PolygonArea \
	AzimuthalEquidistant CassiniSoldner \
	Geoid Gnomonic OSGB AlbersEqualArea
PROGRAMS = GeoConvert TransverseMercatorProj CartConvert Geod GeodesicProj \
	GeoidEval Planimeter ConicProj

HEADERS = $(patsubst %,../include/GeographicLib/%.hpp,Constants $(MODULES))
SOURCES = $(patsubst %,../src/%.cpp,$(MODULES)) \
	$(patsubst %,../tools/%.cpp,$(PROGRAMS))

EXTRAFILES = tmseries30.html geodseries30.html
HTMLMANPAGES = 	$(patsubst %,../man/%.1.html,$(PROGRAMS))
SCRIPTDRIVERS = $(wildcard scripts/*.html)
JSSCRIPTS = $(wildcard scripts/GeographicLib/*.js)

MAXIMA = tm ellint tmseries geod
MAXIMASOURCES = $(patsubst %,../maxima/%.mac,$(MAXIMA))

doc: html/index.html

html/index.html: Doxyfile Geographic.doc \
	$(HEADERS) $(ALLSOURCES) $(MAXIMASOURCES) $(EXTRAFILES) \
	$(HTMLMANPAGES)
	if test -d html; then rm -rf html/*; else mkdir html; fi
	cp -p $(MAXIMASOURCES) $(EXTRAFILES) $(HTMLMANPAGES) html/
	doxygen

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

list:
	@echo Doxyfile Geographic.doc $(EXTRAFILES) \
	$(SCRIPTDRIVERS) $(JSSCRIPTS)

distclean:
	rm -rf html

.PHONY: doc install list clean
