#
# Makefile.am
#
# Copyright (C) 2011, Charles Karney <charles@karney.com>
# $Id$

PYTHON_INTERFACE = GeographicLib.py

PYTHON_ALL = $(PYTHON_INTERFACE)

pythondir=$(datadir)/GeographicLib/python

install:
	$(INSTALL) -d $(DESTDIR)$(pythondir)
	$(INSTALL) -m 644 $(PYTHON_ALL) $(DESTDIR)$(pythondir)

clean-local:
	rm -rf *.pyc

EXTRA_DIST = Makefile.mk CMakeLists.txt $(PYTHON_ALL)
