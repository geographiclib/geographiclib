# $Id: ba2fe05a5521242cf8a105b5686ad7e6e6e11eaa $

PROGRAMS = GeoConvert TransverseMercatorTest CartConvert Geod EquidistantTest \
	GeoidEval Planimeter

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
	@echo GeographicLib-vc8.sln GeographicLib-vc9.sln $(VSPROJECTS)

.PHONY: all install list clean
