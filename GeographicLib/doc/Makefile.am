#  $Id$

EXTRAFILES = tmseries30.html geodseries30.html

HPPFILES = \
$(top_srcdir)/include/GeographicLib/AlbersEqualArea.hpp \
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
$(top_srcdir)/include/GeographicLib/OSGB.hpp \
$(top_srcdir)/include/GeographicLib/PolarStereographic.hpp \
$(top_srcdir)/include/GeographicLib/TransverseMercatorExact.hpp \
$(top_srcdir)/include/GeographicLib/TransverseMercator.hpp \
$(top_srcdir)/include/GeographicLib/UTMUPS.hpp

ALLSOURCES = \
$(top_srcdir)/src/AlbersEqualArea.cpp \
$(top_srcdir)/src/AzimuthalEquidistant.cpp \
$(top_srcdir)/src/CassiniSoldner.cpp \
$(top_srcdir)/src/DMS.cpp \
$(top_srcdir)/src/EllipticFunction.cpp \
$(top_srcdir)/src/Geocentric.cpp \
$(top_srcdir)/src/GeoCoords.cpp \
$(top_srcdir)/src/Geodesic.cpp \
$(top_srcdir)/src/GeodesicLine.cpp \
$(top_srcdir)/src/Geoid.cpp \
$(top_srcdir)/src/Gnomonic.cpp \
$(top_srcdir)/src/LambertConformalConic.cpp \
$(top_srcdir)/src/LocalCartesian.cpp \
$(top_srcdir)/src/MGRS.cpp \
$(top_srcdir)/src/OSGB.cpp \
$(top_srcdir)/src/PolarStereographic.cpp \
$(top_srcdir)/src/TransverseMercator.cpp \
$(top_srcdir)/src/TransverseMercatorExact.cpp \
$(top_srcdir)/src/UTMUPS.cpp \
$(top_srcdir)/tools/CartConvert.cpp \
$(top_srcdir)/tools/EquidistantTest.cpp \
$(top_srcdir)/tools/GeoConvert.cpp \
$(top_srcdir)/tools/Geod.cpp \
$(top_srcdir)/tools/GeoidEval.cpp \
$(top_srcdir)/tools/Planimeter.cpp \
$(top_srcdir)/tools/TransverseMercatorTest.cpp

MANPAGES = \
CartConvert.1.html \
EquidistantTest.1.html \
GeoConvert.1.html \
Geod.1.html \
GeoidEval.1.html \
Planimeter.1.html \
TransverseMercatorTest.1.html

PODFIX= sed -e 's%<head>%<head><link href="http://search.cpan.org/s/style.css" rel="stylesheet" type="text/css">%' -e 's%<code>\([^<>]*\)(\(.\))</code>%<a href="\1.\2.html">&</a>%'g

doc: man html/index.html
man: $(MANPAGES)

CartConvert.1.html: $(top_srcdir)/tools/CartConvert.pod
	pod2html --noindex $^ | $(PODFIX) > $@
EquidistantTest.1.html: $(top_srcdir)/tools/EquidistantTest.pod
	pod2html --noindex $^ | $(PODFIX) > $@
GeoConvert.1.html: $(top_srcdir)/tools/GeoConvert.pod
	pod2html --noindex $^ | $(PODFIX) > $@
Geod.1.html: $(top_srcdir)/tools/Geod.pod
	pod2html --noindex $^ | $(PODFIX) > $@
GeoidEval.1.html: $(top_srcdir)/tools/GeoidEval.pod
	pod2html --noindex $^ | $(PODFIX) > $@
Planimeter.1.html: $(top_srcdir)/tools/Planimeter.pod
	pod2html --noindex $^ | $(PODFIX) > $@
TransverseMercatorTest.1.html: $(top_srcdir)/tools/TransverseMercatorTest.pod
	pod2html --noindex $^ | $(PODFIX) > $@

html/index.html: Doxyfile Geographic.doc \
	$(HPPFILES) $(ALLSOURCES) $(EXTRAFILES) $(MANPAGES)
	if test -d html; then rm -rf html/*; else mkdir html; fi
	cp -p $(EXTRAFILES) $(MANPAGES) $(top_srcdir)/maxima/*.mac html/
	doxygen

distclean-local:
	rm -rf html $(MANPAGES)

install-doc: doc
	$(INSTALL) -d $(DESTDIR)$(docdir)/html
	$(INSTALL) -m 644 html/* $(DESTDIR)$(docdir)/html
