#
# Makefile.am
#
# Copyright (C) 2011, Charles Karney <charles@karney.com>
# $Id$

MANPAGES = \
CartConvert.1 \
ConicProj.1 \
GeodesicProj.1 \
GeoConvert.1 \
Geod.1 \
GeoidEval.1 \
Planimeter.1 \
TransverseMercatorProj.1

USAGE = \
CartConvert.usage \
ConicProj.usage \
GeodesicProj.usage \
GeoConvert.usage \
Geod.usage \
GeoidEval.usage \
Planimeter.usage \
TransverseMercatorProj.usage

HTMLMAN = \
CartConvert.1.html \
ConicProj.1.html \
GeodesicProj.1.html \
GeoConvert.1.html \
Geod.1.html \
GeoidEval.1.html \
Planimeter.1.html \
TransverseMercatorProj.1.html

POD2MAN = pod2man --center="GeographicLib Utilities" \
	--release="$(PACKAGE_STRING)"
PODFIX= sed -e 's%<head>%<head><link href="http://search.cpan.org/s/style.css" rel="stylesheet" type="text/css">%' -e 's%<code>\([^<>]*\)(\(.\))</code>%<a href="\1.\2.html">&</a>%'g
man1_MANS = $(MANPAGES)

SUFFIXES = .pod .1 .usage

all: man
man: manpages usage htmlman

manpages: $(MANPAGES)
usage: $(USAGE)
htmlman: $(HTMLMAN)

.pod.usage:
	sh makeusage.sh $< > $*.usage

.pod.1:
	pod2html --noindex $^ | $(PODFIX) > $*.1.html
	$(POD2MAN) $^ > $@

CartConvert.1.html: CartConvert.pod
	pod2html --noindex $^ | $(PODFIX) > $@
ConicProj.1.html: ConicProj.pod
	pod2html --noindex $^ | $(PODFIX) > $@
GeodesicProj.1.html: GeodesicProj.pod
	pod2html --noindex $^ | $(PODFIX) > $@
GeoConvert.1.html: GeoConvert.pod
	pod2html --noindex $^ | $(PODFIX) > $@
Geod.1.html: Geod.pod
	pod2html --noindex $^ | $(PODFIX) > $@
GeoidEval.1.html: GeoidEval.pod
	pod2html --noindex $^ | $(PODFIX) > $@
Planimeter.1.html: Planimeter.pod
	pod2html --noindex $^ | $(PODFIX) > $@
TransverseMercatorProj.1.html: TransverseMercatorProj.pod
	pod2html --noindex $^ | $(PODFIX) > $@

EXTRA_DIST = Makefile.mk CMakeLists.txt makeusage.sh \
	GeoConvert.pod TransverseMercatorProj.pod \
	CartConvert.pod ConicProj.pod Geod.pod GeodesicProj.pod \
	GeoidEval.pod Planimeter.pod \
	GeoConvert.1 TransverseMercatorProj.1 \
	CartConvert.1 ConicProj.1 Geod.1 GeodesicProj.1 \
	GeoidEval.1 Planimeter.1 \
	GeoConvert.usage TransverseMercatorProj.usage \
	CartConvert.usage ConicProj.usage Geod.usage GeodesicProj.usage \
	GeoidEval.usage Planimeter.usage \
	GeoConvert.1.html TransverseMercatorProj.1.html \
	CartConvert.1.html ConicProj.1.html Geod.1.html GeodesicProj.1.html \
	GeoidEval.1.html Planimeter.1.html

distclean-local:
	rm -rf *.usage *.1.html *.1
