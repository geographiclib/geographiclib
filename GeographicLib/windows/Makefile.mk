# $Id$

PROGRAMS = GeoConvert TransverseMercatorTest CartConvert Geod EquidistantTest \
	GeoidEval Planimeter

VSPROJECTS = $(addsuffix -vc8.vcproj,Geographic $(PROGRAMS)) \
	$(addsuffix -vc9.vcproj,Geographic $(PROGRAMS))

all:
	@:
install:
	@:
clean:
	@:
list:
	@echo GeographicLib-vc8.sln GeographicLib-vc9.sln $(VSPROJECTS)

.PHONY: all install list clean
