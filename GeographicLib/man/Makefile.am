#
# Makefile.am
#
# Copyright (C) 2011, Charles Karney <charles@karney.com>


MANPAGES = \
CartConvert.1 \
EquidistantTest.1 \
GeoConvert.1 \
Geod.1 \
GeoidEval.1 \
Planimeter.1 \
TransverseMercatorTest.1

POD2MAN = pod2man --center="GeographicLib Utilities"
man1_MANS = $(MANPAGES)

CartConvert.1: $(top_srcdir)/tools/CartConvert.pod
	$(POD2MAN) $^ > $@
EquidistantTest.1: $(top_srcdir)/tools/EquidistantTest.pod
	$(POD2MAN) $^ > $@
GeoConvert.1: $(top_srcdir)/tools/GeoConvert.pod
	$(POD2MAN) $^ > $@
Geod.1: $(top_srcdir)/tools/Geod.pod
	$(POD2MAN) $^ > $@
GeoidEval.1: $(top_srcdir)/tools/GeoidEval.pod
	$(POD2MAN) $^ > $@
Planimeter.1: $(top_srcdir)/tools/Planimeter.pod
	$(POD2MAN) $^ > $@
TransverseMercatorTest.1: $(top_srcdir)/tools/TransverseMercatorTest.pod
	$(POD2MAN) $^ > $@
