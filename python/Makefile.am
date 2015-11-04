#
# Makefile.am
#
# Copyright (C) 2011, Charles Karney <charles@karney.com>

PACKAGE=geographiclib
PYTHON_FILES = \
	$(srcdir)/$(PACKAGE)/__init__.py \
	$(srcdir)/$(PACKAGE)/geomath.py \
	$(srcdir)/$(PACKAGE)/constants.py \
	$(srcdir)/$(PACKAGE)/accumulator.py \
	$(srcdir)/$(PACKAGE)/geodesiccapability.py \
	$(srcdir)/$(PACKAGE)/geodesic.py \
	$(srcdir)/$(PACKAGE)/geodesicline.py \
	$(srcdir)/$(PACKAGE)/polygonarea.py
TEST_FILES = \
	$(srcdir)/$(PACKAGE)/test/__init__.py \
	$(srcdir)/$(PACKAGE)/test/test_geodesic.py

pythondir=$(libdir)/python/site-packages/$(PACKAGE)

install:
	$(INSTALL) -d $(DESTDIR)$(pythondir)
	$(INSTALL) -m 644 $(PYTHON_FILES) $(DESTDIR)$(pythondir)
	$(INSTALL) -d $(DESTDIR)$(pythondir)/test
	$(INSTALL) -m 644 $(TEST_FILES) $(DESTDIR)$(pythondir)/test

# Don't install setup.py because it ends up in e.g.,
# /usr/local/lib/python/site-packages/setup.py
#	$(INSTALL) -m 644 setup.py $(DESTDIR)$(pythondir)/../

clean-local:
	rm -rf *.pyc $(PACKAGE)/*.pyc

EXTRA_DIST = Makefile.mk $(PACKAGE)/CMakeLists.txt $(PYTHON_FILES) setup.py \
	MANIFEST.in README.txt
