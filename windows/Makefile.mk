# $Id$

PROGRAMS = GeoConvert \
	TransverseMercatorProj \
	CartConvert \
	Geod \
	GeodesicProj \
	GeoidEval \
	Planimeter \
	ConicProj

VSPROJECTS = $(addsuffix -vc8.vcproj,Geographic $(PROGRAMS)) \
	$(addsuffix -vc9.vcproj,Geographic $(PROGRAMS)) \
	$(addsuffix -vc10.vcxproj,Geographic $(PROGRAMS))

all:
	@:
install:
	@:
clean:
	@:

.PHONY: all install clean
