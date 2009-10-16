# $Id$

MAKE := $(MAKE) -f $(lastword $(MAKEFILE_LIST))
SUBDIRS = src tools doc
ALLDIRS = include $(SUBDIRS) maxima windows

all: src tools

$(SUBDIRS):
	$(MAKE) -C $@

tools: src
install: install-headers install-lib install-tools
clean: clean-src clean-tools clean-doc
install-headers:
	$(MAKE) -C include install
install-lib:
	$(MAKE) -C src install
install-tools: src
	$(MAKE) -C tools install
install-doc: doc
	$(MAKE) -C doc install
clean-src:
	$(MAKE) -C src clean
clean-tools:
	$(MAKE) -C tools clean
clean-doc:
	$(MAKE) -C doc clean

list:
	@echo 00README.txt; echo COPYING; echo AUTHORS; echo Makefile
	@for d in $(ALLDIRS); do \
	  (echo Makefile; $(MAKE) -s -C $$d list) | tr ' ' '\n' | \
	  while read f; do echo $$d/$$f; done; \
	done

package:
	test -d distrib || mkdir distrib
	$(MAKE) -s list | while read f;do \
	  echo GeographicLib/$$f; \
	done | xargs tar Ccfzv .. distrib/Geographic.tgz

.PHONY: all $(SUBDIRS) \
	install install-headers install-lib install-tools install-doc \
	clean clean-src clean-tools clean-doc list package
