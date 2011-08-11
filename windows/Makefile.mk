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
list:
	@echo $(VSPROJECTS) \
	GeographicLib-vc8.sln GeographicLib-vc9.sln GeographicLib-vc10.sln 

.PHONY: all install list clean
