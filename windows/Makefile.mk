# $Id: 959eb5fb059ac028bdef404e4cb187a4a81bb318 $

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
