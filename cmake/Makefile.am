#
# Makefile.am
#
# Copyright (C) 2011, Charles Karney <charles@karney.com>
# $Id: a072d10e7ea813c07a37d1d612a72ab82d0406d2 $

cmakedir=$(datadir)/cmake/GeographicLib

install:
	$(INSTALL) -d $(DESTDIR)$(cmakedir)
	$(INSTALL) -m 644 FindGeographicLib.cmake $(DESTDIR)$(cmakedir)

EXTRA_DIST = Makefile.mk FindGeographicLib.cmake \
	geographiclib-config-version.cmake.in geographiclib-config.cmake.in
