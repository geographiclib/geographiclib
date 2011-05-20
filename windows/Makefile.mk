# $Id: Makefile.mk 6768 2009-12-04 15:39:30Z karney $

PROGRAMS = GeoConvert TransverseMercatorTest CartConvert Geod EquidistantTest \
	GeoidEval

VSPROJECTS = $(addsuffix -vc8.vcproj,GeographicLib $(PROGRAMS)) \
	$(addsuffix -vc9.vcproj,GeographicLib $(PROGRAMS))

all:
	@:
install:
	@:
clean:
	@:
list:
	@echo GeographicLib-vc8.sln GeographicLib-vc9.sln $(VSPROJECTS)

.PHONY: all install list clean
