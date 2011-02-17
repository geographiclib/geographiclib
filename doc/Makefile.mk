# $Id$

MODULES = DMS EllipticFunction GeoCoords MGRS PolarStereographic \
	TransverseMercator TransverseMercatorExact UTMUPS Geocentric \
	LocalCartesian Geodesic GeodesicLine \
	AzimuthalEquidistant CassiniSoldner \
	Geoid Gnomonic OSGB AlbersEqualArea
PROGRAMS = GeoConvert TransverseMercatorTest CartConvert Geod EquidistantTest \
	GeoidEval Planimeter

HEADERS = Constants.hpp $(patsubst %,../include/GeographicLib/%.hpp,$(MODULES))
SOURCES = $(patsubst %,../src/%.cpp,$(MODULES)) \
	$(patsubst %,../tools/%.cpp,$(PROGRAMS))

EXTRAFILES = tmseries30.html geodseries30.html
HTMLMANPAGES = $(addsuffix .1.html,$(PROGRAMS))

%.1.html: %.pod
	pod2html --noindex $^ | sed -e 's%<head>%<head><link href="http://search.cpan.org/s/style.css" rel="stylesheet" type="text/css">%' -e 's%<code>\([^<>]*\)(\(.\))</code>%<a href="\1.\2.html">&</a>%'g > $@

MAXIMA = tm ellint tmseries geod
MAXIMASOURCES = $(patsubst %,../maxima/%.mac,$(MAXIMA))

doc: man html/index.html

VPATH = ../src ../include/GeographicLib ../tools ../maxima

man: $(HTMLMANPAGES)

html/index.html: Doxyfile Geographic.doc \
	$(HEADERS) $(ALLSOURCES) $(MAXIMASOURCES) $(EXTRAFILES) \
	$(HTMLMANPAGES)
	if test -d html; then rm -rf html/*; else mkdir html; fi
	cp -p $(MAXIMAASOURCES) $(EXTRAFILES) $(HTMLMANPAGES) html/
	doxygen

PREFIX = /usr/local
DEST = $(PREFIX)/share/GeographicLib/doc/html
INSTALL = install -b

install: man html/index.html
	test -d $(DEST) || mkdir -p $(DEST)
	$(INSTALL) -m 644 html/* $(DEST)/
list:
	@echo Doxyfile Geographic.doc $(EXTRAFILES)
clean:
	rm -rf html

.PHONY: doc install list clean
