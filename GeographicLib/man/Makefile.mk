# $Id$
PROGRAMS = GeoConvert TransverseMercatorTest CartConvert Geod EquidistantTest \
	GeoidEval Planimeter

PREFIX = /usr/local
MANPAGES = $(addsuffix .1,$(PROGRAMS))
DEST = $(PREFIX)/share/man/man1

VERSION:=$(shell grep '\bVERSION=' ../configure | cut -f2 -d\' | head -1)

VPATH = ../tools

%.1: %.pod
	pod2man --center="GeographicLib Utilities" \
	--release="GeographicLib $(VERSION)" $^ > $@

all: $(MANPAGES)

INSTALL=install -b

install:
	test -d $(DEST) || mkdir -p $(DEST)
	$(INSTALL) -m 644 $(MANPAGES) $(DEST)/
list:
	@echo $(MANPAGES)

distclean:
	rm -f *.1

.PHONY: all install list clean
