# $Id$

MAKEFILE := $(lastword $(MAKEFILE_LIST))
MAKE := $(MAKE) -f $(MAKEFILE)
SUBDIRS = src man tools doc
ALLDIRS = include $(SUBDIRS) maxima matlab windows

all: src man tools

$(SUBDIRS):
	$(MAKE) -C $@

tools: src
install: install-headers install-lib install-tools install-man
clean: clean-src clean-tools clean-doc clean-man
distclean: clean distclean-doc distclean-man
install-headers:
	$(MAKE) -C include install
install-lib:
	$(MAKE) -C src install
install-tools: src
	$(MAKE) -C tools install
install-doc: doc
	$(MAKE) -C doc install
install-man: man
	$(MAKE) -C man install
clean-src:
	$(MAKE) -C src clean
clean-tools:
	$(MAKE) -C tools clean
clean-doc:
	$(MAKE) -C doc clean
clean-man:
	$(MAKE) -C man clean
distclean-doc:
	$(MAKE) -C doc distclean
distclean-man:
	$(MAKE) -C man distclean

list:
	@for f in 00README.txt COPYING.txt AUTHORS NEWS Makefile \
	$(MAKEFILE); do \
	  echo $$f; \
	done
	@for d in $(ALLDIRS); do \
	  (echo $(MAKEFILE); $(MAKE) -s -C $$d list) | tr ' ' '\n' | \
	  while read f; do echo $$d/$$f; done; \
	done

VERSION:=$(shell grep '\bVERSION=' configure | cut -f2 -d\' | head -1)

.PHONY: all $(SUBDIRS) \
	install install-headers install-lib install-tools install-doc \
	clean clean-src clean-tools clean-doc list package
