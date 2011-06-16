# $Id$

PROGRAMS = GeoConvert \
	TransverseMercatorProj \
	CartConvert \
	Geod \
	GeodesicProj \
	GeoidEval \
	Planimeter \
	ConicProj
SCRIPTS = geographiclib-get-geoids

all: $(PROGRAMS) $(SCRIPTS)

LIBSTEM = Geographic
LIBRARY = lib$(LIBSTEM).a

INCLUDEPATH = ../include
LIBPATH = ../src

# After installation, use these values of INCLUDEPATH and LIBPATH
# INCLUDEPATH = $(PREFIX)/include
# LIBPATH = $(PREFIX)/lib

PREFIX = /usr/local
GEOID_DEFAULT_PATH = $(PREFIX)/share/GeographicLib/geoids

CC = g++ -g
CXXFLAGS = -g -Wall -Wextra -O3

CPPFLAGS = -I$(INCLUDEPATH) -I../man $(DEFINES)
LDLIBS = -L$(LIBPATH) -l$(LIBSTEM)

$(PROGRAMS): $(LIBPATH)/$(LIBRARY)
	$(CC) $(LDFLAGS) -o $@ $@.o $(LDLIBS)

VPATH = ../include/GeographicLib ../man

clean:
	rm -f *.o $(SCRIPTS)

GeoConvert: GeoConvert.o
TransverseMercatorProj: TransverseMercatorProj.o
CartConvert: CartConvert.o
Geod: Geod.o
GeodesicProj: GeodesicProj.o
GeoidEval: GeoidEval.o
Planimeter: Planimeter.o
ConicProj: ConicProj.o

GeoConvert.o: GeoConvert.usage Constants.hpp Config.h DMS.hpp \
	GeoCoords.hpp UTMUPS.hpp
TransverseMercatorProj.o: TransverseMercatorProj.usage Constants.hpp \
	Config.h DMS.hpp EllipticFunction.hpp TransverseMercator.hpp \
	TransverseMercatorExact.hpp
CartConvert.o: CartConvert.usage Constants.hpp Config.h DMS.hpp \
	Geocentric.hpp LocalCartesian.hpp
Geod.o: Geod.usage Constants.hpp Config.h DMS.hpp Geodesic.hpp \
	GeodesicLine.hpp
GeodesicProj.o: GeodesicProj.usage AzimuthalEquidistant.hpp \
	CassiniSoldner.hpp Gnomonic.hpp Constants.hpp Config.h DMS.hpp \
	Geodesic.hpp GeodesicLine.hpp
GeoidEval.o: GeoidEval.usage Constants.hpp Config.h DMS.hpp \
	GeoCoords.hpp Geoid.hpp
Planimeter.o: Planimeter.usage Constants.hpp Config.h DMS.hpp \
	GeoCoords.hpp Geodesic.hpp GeodesicLine.hpp
ConicProj.o: ConicProj.usage AlbersEqualArea.hpp \
	Constants.hpp Config.h DMS.hpp LambertConformalConic.hpp

geographiclib-get-geoids: geographiclib-get-geoids.sh
	sed -e "s%@GEOID_DEFAULT_PATH@%$(GEOID_DEFAULT_PATH)%" $< > $@
	chmod +x $@

INSTALL = install -b
DEST = $(PREFIX)/bin
SDEST = $(PREFIX)/sbin

list:
	@echo $(addsuffix .cpp,$(PROGRAMS))

install: $(PROGRAMS) $(SCRIPTS)
	test -f $(DEST) || mkdir -p $(DEST)
	$(INSTALL) $(PROGRAMS) $(DEST)
	test -f $(SDEST) || mkdir -p $(SDEST)
	$(INSTALL) $(SCRIPTS) $(SDEST)
