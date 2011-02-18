# $Id$
PROGRAMS = GeoConvert TransverseMercatorTest CartConvert Geod EquidistantTest \
	GeoidEval Planimeter

PREFIX = /usr/local
MANPAGES = $(addsuffix .1,$(PROGRAMS))
DEST = $(PREFIX)/share/man/man1

VPATH = ../tools

%.1: %.pod
	pod2man --center="GeographicLib Utilities" $^ > $@

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
