#  $Id: Makefile.am 6859 2010-09-06 14:45:33Z karney $

FIGURES = gauss-krueger-graticule thompson-tm-graticule \
	gauss-krueger-convergence-scale gauss-schreiber-graticule-a \
	gauss-krueger-graticule-a thompson-tm-graticule-a \
	gauss-krueger-error

EXTRAFILES = tmseries30.html geodseries30.html

HPPFILES = \
$(top_srcdir)/include/GeographicLib/AzimuthalEquidistant.hpp \
$(top_srcdir)/include/GeographicLib/CassiniSoldner.hpp \
$(top_srcdir)/include/GeographicLib/Constants.hpp \
$(top_srcdir)/include/GeographicLib/DMS.hpp \
$(top_srcdir)/include/GeographicLib/EllipticFunction.hpp \
$(top_srcdir)/include/GeographicLib/Geocentric.hpp \
$(top_srcdir)/include/GeographicLib/GeoCoords.hpp \
$(top_srcdir)/include/GeographicLib/Geodesic.hpp \
$(top_srcdir)/include/GeographicLib/GeodesicLine.hpp \
$(top_srcdir)/include/GeographicLib/Geoid.hpp \
$(top_srcdir)/include/GeographicLib/Gnomonic.hpp \
$(top_srcdir)/include/GeographicLib/LambertConformalConic.hpp \
$(top_srcdir)/include/GeographicLib/LocalCartesian.hpp \
$(top_srcdir)/include/GeographicLib/MGRS.hpp \
$(top_srcdir)/include/GeographicLib/PolarStereographic.hpp \
$(top_srcdir)/include/GeographicLib/TransverseMercatorExact.hpp \
$(top_srcdir)/include/GeographicLib/TransverseMercator.hpp \
$(top_srcdir)/include/GeographicLib/UTMUPS.hpp

ALLSOURCES = \
$(top_srcdir)/src/AzimuthalEquidistant.cpp \
$(top_srcdir)/src/CassiniSoldner.cpp \
$(top_srcdir)/src/DMS.cpp \
$(top_srcdir)/src/EllipticFunction.cpp \
$(top_srcdir)/src/Geocentric.cpp \
$(top_srcdir)/src/GeoCoords.cpp \
$(top_srcdir)/src/Geodesic.cpp \
$(top_srcdir)/src/Geoid.cpp \
$(top_srcdir)/src/Gnomonic.cpp \
$(top_srcdir)/src/LambertConformalConic.cpp \
$(top_srcdir)/src/LocalCartesian.cpp \
$(top_srcdir)/src/MGRS.cpp \
$(top_srcdir)/src/PolarStereographic.cpp \
$(top_srcdir)/src/TransverseMercator.cpp \
$(top_srcdir)/src/TransverseMercatorExact.cpp \
$(top_srcdir)/src/UTMUPS.cpp \
$(top_srcdir)/tools/CartConvert.cpp \
$(top_srcdir)/tools/EquidistantTest.cpp \
$(top_srcdir)/tools/GeoConvert.cpp \
$(top_srcdir)/tools/Geod.cpp \
$(top_srcdir)/tools/GeoidEval.cpp \
$(top_srcdir)/tools/TransverseMercatorTest.cpp

doc: html/index.html
html/index.html: Doxyfile Geographic.doc \
	$(HPPFILES) $(ALLSOURCES) \
	$(addsuffix .pdf,$(FIGURES)) $(addsuffix .png,$(FIGURES)) \
	$(EXTRAFILES)
	rm -rf html/*
	doxygen
	for f in $(FIGURES); do cp -p $$f.pdf $$f.png html/; done
	for f in $(EXTRAFILES); do cp -p $$f html/; done
	for f in $(top_srcdir)/maxima/*.mac; do cp -p $$f html/.; done

clean-local:
	rm -rf html

install-doc: doc
	$(INSTALL) -d $(DESTDIR)$(docdir)/html
	$(INSTALL) -m 644 html/* $(DESTDIR)$(docdir)/html

