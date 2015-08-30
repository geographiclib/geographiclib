JS_MODULES=Math Geodesic GeodesicLine PolygonArea DMS Interface
JS_SCRIPTS = $(patsubst %,GeographicLib/%.js,$(JS_MODULES))

SCRIPTDRIVERSIN = $(wildcard geod-*.in)
SCRIPTDRIVERS = $(patsubst %.in,%.html,$(SCRIPTDRIVERSIN))

all: geographiclib.js $(SCRIPTDRIVERS)

%.html: %.in
	cp $^ $@

geographiclib.js: HEADER.js $(JS_SCRIPTS)
	./js-compress.sh $^ > $@

clean:
	rm -f geographiclib.js *.html

PREFIX = /usr/local
DEST = $(PREFIX)/share/doc/GeographicLib/scripts
INSTALL = install -b

install: all
	test -d $(DEST)/GeographicLib || \
	mkdir -p $(DEST)/GeographicLib
	$(INSTALL) -m 644 $(SCRIPTDRIVERS) geographiclib.js $(DEST)/
	$(INSTALL) -m 644 $(JS_SCRIPTS) $(DEST)/GeographicLib/

.PHONY: install clean
