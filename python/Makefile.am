#
# Makefile.am
#
# Copyright (C) 2011, Charles Karney <charles@karney.com>
# $Id$

PACKAGE=geographiclib
PYTHON_FILES = \
	$(PACKAGE)/__init__.py \
	$(PACKAGE)/geomath.py \
	$(PACKAGE)/constants.py \
	$(PACKAGE)/accumulator.py \
	$(PACKAGE)/geodesiccapability.py \
	$(PACKAGE)/geodesic.py \
	$(PACKAGE)/geodesicline.py \
	$(PACKAGE)/polygonarea.py

pythondir=$(libdir)/python/site-packages/$(PACKAGE)

install:
	$(INSTALL) -d $(DESTDIR)$(pythondir)
	$(INSTALL) -m 644 $(PYTHON_FILES) $(DESTDIR)$(pythondir)
	$(INSTALL) -m 644 setup.py $(DESTDIR)$(pythondir)/../

clean-local:
	rm -rf $(PACKAGE)/*.pyc

EXTRA_DIST = Makefile.mk $(PACKAGE)/CMakeLists.txt $(PYTHON_FILES) setup.py
