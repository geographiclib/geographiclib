# $Id: Makefile.mk 6722 2009-10-18 15:28:50Z ckarney $

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
