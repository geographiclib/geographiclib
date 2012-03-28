#
# Makefile.am
#
# Copyright (C) 2011, Charles Karney <charles@karney.com>
# $Id$

cmakedir=$(datadir)/cmake/GeographicLib

install:
	$(INSTALL) -d $(DESTDIR)$(cmakedir)
	$(INSTALL) -m 644 $(srcdir)/FindGeographicLib.cmake \
		$(DESTDIR)$(cmakedir)

EXTRA_DIST = Makefile.mk FindGeographicLib.cmake \
	geographiclib-config-version.cmake.in geographiclib-config.cmake.in
