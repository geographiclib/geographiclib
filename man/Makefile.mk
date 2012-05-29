PROGRAMS = CartConvert \
	ConicProj \
	GeoConvert \
	Geod \
	GeodesicProj \
	GeoidEval \
	Gravity \
	MagneticField \
	Planimeter \
	TransverseMercatorProj

PODSRC = $(addsuffix .pod,$(PROGRAMS))
MANPAGES = $(addsuffix .1,$(PROGRAMS))
USAGE = $(addsuffix .usage,$(PROGRAMS))
HTMLMAN = $(addsuffix .1.html,$(PROGRAMS))

PREFIX = /usr/local

DEST = $(PREFIX)/share/man/man1

VERSION:=$(shell grep '\bVERSION=' ../configure | cut -f2 -d\' | head -1)

%.1: %.pod
	pod2man --center="GeographicLib Utilities" \
	--release="GeographicLib $(VERSION)" $^ > $@

%.1.html: %.pod
	pod2html --noindex $^ | sed -e 's%<head>%<head><link href="http://search.cpan.org/s/style.css" rel="stylesheet" type="text/css">%' -e 's%<code>\([^<>]*\)(\(.\))</code>%<a href="\1.\2.html">&</a>%'g > $@

%.usage: %.pod
	sh makeusage.sh $< > $@

all: $(MANPAGES) $(USAGE) $(HTMLMAN)

INSTALL = install -b

install:
	test -d $(DEST) || mkdir -p $(DEST)
	$(INSTALL) -m 644 $(MANPAGES) $(DEST)/

maintainer-clean:
	rm -f *.1 *.usage *.1.html

.PHONY: all install clean
