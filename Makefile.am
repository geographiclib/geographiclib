#
# Makefile.am
#
# Copyright (C) 2009, Francesco P. Lovergine <frankie@debian.org>
# $Id$

AUTOMAKE_OPTIONS = foreign

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src man tools doc include matlab python

EXTRA_DIST = AUTHORS 00README.txt LICENSE.txt NEWS INSTALL \
	Makefile.mk CMakeLists.txt \
	cmake windows maxima man doc python

dist-hook:
	find $(distdir)/maxima -type f -name '*.lsp' | xargs rm -rf ; \
	find $(distdir)/cmake -type f \
		! \(-name '*.cmake' -o -name '*.mk' \) | xargs rm -rf ; \
	find $(distdir)/windows -mindepth 1 -type d | xargs rm -rf ; \
	find $(distdir)/windows -type f \
		! \( -name '*.sln' -o -name '*.vc*proj' -o -name '*.mk' \)| \
		xargs rm -f ; \
	find $(distdir)/matlab -type f \
		! \( -name '*.cpp' -o -name '*.m' -o -name 'Makefile.*' -o -name CMakeLists.txt \) | \
		xargs rm -rf ; \
	find $(distdir) \
		\( -name .svn -o -name '.git*' -o -name CVS -o -name Makefile -o -name '*~' -o -name '*#*' -o -name 'CMakeFiles' -o -name '*.log' -o -name '*.cmake' -o -name '*.tmp' -o -name '*.pyc' \)| \
		xargs rm -rf ; \
	echo include Makefile.mk > $(distdir)/Makefile
	echo '#define' GEOGRAPHICLIB_VERSION_STRING \
	\"$(PACKAGE_VERSION)\" > $(distdir)/include/GeographicLib/Config.h
# Custom rules

all-local: man doc
install-data-local: install-doc # install-matlab

doc: man
	$(MAKE) -C doc doc

install-doc:
	$(MAKE) -C doc install-doc

man:
	$(MAKE) -C man man

# install-matlab:
# 	$(MAKE) -C matlab install-matlab

.PHONY: doc install-doc man install-matlab install-python
