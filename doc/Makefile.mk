# $Id: Makefile.mk 6859 2010-09-06 14:45:33Z karney $

MODULES = DMS EllipticFunction GeoCoords MGRS PolarStereographic \
	TransverseMercator TransverseMercatorExact UTMUPS Geocentric \
	LocalCartesian Geodesic GeodesicLine \
	AzimuthalEquidistant CassiniSoldner \
	Geoid Gnomonic
PROGRAMS = GeoConvert TransverseMercatorTest CartConvert Geod EquidistantTest \
	GeoidEval

HEADERS = Constants.hpp $(patsubst %,../include/GeographicLib/%.hpp,$(MODULES))
SOURCES = $(patsubst %,../src/%.cpp,$(MODULES)) \
	$(patsubst %,../tools/%.cpp,$(PROGRAMS))

FIGURES = gauss-krueger-graticule thompson-tm-graticule \
	gauss-krueger-convergence-scale gauss-schreiber-graticule-a \
	gauss-krueger-graticule-a thompson-tm-graticule-a \
	gauss-krueger-error
FIGURESOURCES = $(addsuffix .pdf,$(FIGURES)) $(addsuffix .png,$(FIGURES))

EXTRAFILES = tmseries30.html geodseries30.html

MAXIMA = tm ellint tmseries geod
MAXIMASOURCES = $(patsubst %,../maxima/%.mac,$(MAXIMA))

doc: html/index.html

VPATH = ../src ../include/GeographicLib ../tools ../maxima

html/index.html: Doxyfile Geographic.doc \
	$(HEADERS) $(ALLSOURCES) $(FIGURESOURCES) $(MAXIMASOURCES) $(EXTRAFILES)
	if test -d html; then rm -rf html/*; else mkdir html; fi
	for f in $(FIGURESOURCES); do cp -p $$f html/; done
	for f in $(MAXIMASOURCES); do cp -p $$f html/; done
	for f in $(EXTRAFILES); do cp -p $$f html/; done
	doxygen

PREFIX = /usr/local
DEST = $(PREFIX)/share/GeographicLib/doc/html
INSTALL = install -b

install: html/index.html
	test -d $(DEST) || mkdir -p $(DEST)
	$(INSTALL) -m 644 html/* $(DEST)/
list:
	@echo Doxyfile Geographic.doc $(FIGURESOURCES) $(EXTRAFILES)
clean:
	rm -rf html

.PHONY: doc install list clean
