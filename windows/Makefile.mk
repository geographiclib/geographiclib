# $Id: Makefile.mk 6869 2010-09-22 15:50:45Z karney $

PROGRAMS = GeoConvert TransverseMercatorTest CartConvert Geod EquidistantTest \
	GeoidEval Planimeter

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
