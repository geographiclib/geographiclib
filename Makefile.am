#
# Makefile.am
#
# Copyright (C) 2009, Francesco P. Lovergine <frankie@debian.org>
# $Id: Makefile.am 6845 2010-07-21 10:23:46Z karney $

AUTOMAKE_OPTIONS = foreign

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src tools doc include

EXTRA_DIST = AUTHORS 00README.txt COPYING NEWS INSTALL \
			 Makefile.mk windows maxima doc 

dist-hook:
	find $(distdir)/maxima -type f -name '*.lsp' | xargs rm -rf ; \
	find $(distdir)/windows -mindepth 1 -type d | xargs rm -rf ; \
	find $(distdir)/windows -type f \
		\( ! -name '*.sln' -a ! -name '*.vcproj' -a ! -name '*.mk' \)| \
		xargs rm -f ; \
	find $(distdir) \
		\( -name .svn -o -name .git -o -name CVS -o -name Makefile -o -name '*~' \)| \
		xargs rm -rf ; \
	echo include Makefile.mk > $(distdir)/Makefile

# Custom rules

all-local: doc
install-data-local: install-doc

doc:
	$(MAKE) -C doc doc

install-doc:
	$(MAKE) -C doc install-doc

.PHONY: doc install-doc
