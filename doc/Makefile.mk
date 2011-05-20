# $Id$

MODULES = DMS EllipticFunction GeoCoords MGRS PolarStereographic \
	TransverseMercator TransverseMercatorExact UTMUPS Geocentric \
	LocalCartesian Geodesic GeodesicLine \
	AzimuthalEquidistant CassiniSoldner \
	Geoid Gnomonic OSGB AlbersEqualArea
PROGRAMS = GeoConvert TransverseMercatorProj CartConvert Geod GeodesicProj \
	GeoidEval Planimeter ConicProj

HEADERS = $(patsubst %,../include/GeographicLib/%.hpp,Constants $(MODULES))
SOURCES = $(patsubst %,../src/%.cpp,$(MODULES)) \
	$(patsubst %,../tools/%.cpp,$(PROGRAMS))

EXTRAFILES = tmseries30.html geodseries30.html
HTMLMANPAGES = 	$(patsubst %,../man/%.1.html,$(PROGRAMS))

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
DEST = $(PREFIX)/share/doc/GeographicLib/html
INSTALL = install -b

install: html/index.html
	test -d $(DEST) || mkdir -p $(DEST)
	$(INSTALL) -m 644 html/* $(DEST)/
list:
	@echo Doxyfile Geographic.doc $(EXTRAFILES)

distclean:
	rm -rf html

.PHONY: doc install list clean
