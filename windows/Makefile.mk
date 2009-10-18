# $Id$

PROGRAMS = GeoConvert TransverseMercatorTest CartConvert Geod EquidistantTest \
	GeoidEval

VSPROJECTS = $(addsuffix .vcproj,GeographicLib $(PROGRAMS))

all:
	@:
install:
	@:
clean:
	@:
list:
	@echo GeographicLib.sln $(VSPROJECTS)

.PHONY: all install list clean
