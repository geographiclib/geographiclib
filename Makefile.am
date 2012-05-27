#
# Makefile.am
#
# Copyright (C) 2009, Francesco P. Lovergine <frankie@debian.org>
# $Id$

AUTOMAKE_OPTIONS = foreign

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src man tools doc include matlab python cmake examples

EXTRA_DIST = AUTHORS 00README.txt LICENSE.txt NEWS INSTALL \
	Makefile.mk CMakeLists.txt windows maxima doc

dist-hook:
	find $(distdir)/maxima -type f -name '*.lsp' | xargs rm -rf ; \
	find $(distdir)/windows -mindepth 1 -type d | xargs rm -rf ; \
	find $(distdir)/windows -type f \
		! \( -name '*.sln' -o -name '*.vc*proj' -o -name '*.mk' \)| \
		xargs rm -f ; \
	find $(distdir) \
		\( -name .svn -o -name '.git*' -o -name CVS -o -name Makefile -o -name '*~' -o -name '*#*' -o -name 'CMakeFiles' -o -name '*.log' -o -name '*.tmp' -o -name '*.pyc' -o -name '*.bak' -o -name '*.BAK' -o -name geographiclib.js \)| \
		xargs rm -rf ; \
	echo include Makefile.mk > $(distdir)/Makefile ; \
	sed -e "s/Unconfigured/$(PACKAGE_VERSION)/" \
		include/GeographicLib/Config.h > \
		$(distdir)/include/GeographicLib/Config.h
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
