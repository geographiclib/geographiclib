# $Id: Makefile.mk 6816 2010-02-05 21:03:10Z karney $

MAKEFILE := $(lastword $(MAKEFILE_LIST))
MAKE := $(MAKE) -f $(MAKEFILE)
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
	@for f in 00README.txt COPYING AUTHORS NEWS Makefile $(MAKEFILE); do \
	  echo $$f; \
	done
	@for d in $(ALLDIRS); do \
	  (echo $(MAKEFILE); $(MAKE) -s -C $$d list) | tr ' ' '\n' | \
	  while read f; do echo $$d/$$f; done; \
	done

VERSION:=$(shell grep '\bVERSION=' configure | cut -f2 -d\')

package:
	echo include Makefile.mk > Makefile
	test -d distrib || mkdir distrib
	$(MAKE) -s list | while read f;do \
	  echo GeographicLib/$$f; \
	done | xargs tar Ccfz .. distrib/Geographic-$(VERSION).tgz
	rm -rf distrib/GeographicLib
	tar Cxfpz distrib distrib/Geographic-$(VERSION).tgz
	cd distrib && zip -r Geographic-$(VERSION).zip GeographicLib && \
	  rm -rf GeographicLib

.PHONY: all $(SUBDIRS) \
	install install-headers install-lib install-tools install-doc \
	clean clean-src clean-tools clean-doc list package
