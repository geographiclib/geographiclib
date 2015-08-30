EXTRAFILES = $(srcdir)/HEADER.js

SCRIPTDRIVERS = \
	geod-calc.html \
	geod-google.html \
	geod-google-instructions.html

JSSCRIPTS = \
	$(srcdir)/GeographicLib/Math.js \
	$(srcdir)/GeographicLib/Geodesic.js \
	$(srcdir)/GeographicLib/GeodesicLine.js \
	$(srcdir)/GeographicLib/PolygonArea.js \
	$(srcdir)/GeographicLib/DMS.js \
	$(srcdir)/GeographicLib/Interface.js

all: geographiclib.js $(SCRIPTDRIVERS)

.in.html:
	cp $^ $@

geographiclib.js: HEADER.js $(JSSCRIPTS)
	$(srcdir)/js-compress.sh $^ > $@

maintainer-clean-local:
	rm -rf html manpages

scriptsdir=$(DESTDIR)$(docdir)/scripts

install:
	$(INSTALL) -d $(scriptsdir)
	$(INSTALL) -m 644 $(SCRIPTDRIVERS) geographiclib.js $(scriptsdir)
	$(INSTALL) -d $(scriptsdir)/GeographicLib
	$(INSTALL) -m 644 $(JSSCRIPTS) $(scriptsdir)/GeographicLib
