# $Id$

PROGRAMS = GeoConvert \
	TransverseMercatorProj \
	CartConvert \
	Geod \
	GeodesicProj \
	GeoidEval \
	MagneticField \
	Planimeter \
	ConicProj
SCRIPTS = geographiclib-get-geoids \
	geographiclib-get-gravity \
	geographiclib-get-magnetic

all: $(PROGRAMS) $(SCRIPTS)

LIBSTEM = Geographic
LIBRARY = lib$(LIBSTEM).a

INCLUDEPATH = ../include
LIBPATH = ../src

# After installation, use these values of INCLUDEPATH and LIBPATH
# INCLUDEPATH = $(PREFIX)/include
# LIBPATH = $(PREFIX)/lib

PREFIX = /usr/local
GEOGRAPHICLIB_DATA = $(PREFIX)/share/GeographicLib

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
MagneticField: MagneticField.o
Planimeter: Planimeter.o
ConicProj: ConicProj.o

CartConvert.o: CartConvert.usage Config.h Constants.hpp DMS.hpp Geocentric.hpp \
	LocalCartesian.hpp Math.hpp Utility.hpp
ConicProj.o: ConicProj.usage AlbersEqualArea.hpp Config.h Constants.hpp \
	DMS.hpp LambertConformalConic.hpp Math.hpp Utility.hpp
GeoConvert.o: GeoConvert.usage Config.h Constants.hpp DMS.hpp GeoCoords.hpp \
	Math.hpp UTMUPS.hpp Utility.hpp
Geod.o: Geod.usage Config.h Constants.hpp DMS.hpp Geodesic.hpp \
	GeodesicLine.hpp Math.hpp Utility.hpp
GeodesicProj.o: GeodesicProj.usage AzimuthalEquidistant.hpp CassiniSoldner.hpp \
	Config.h Constants.hpp DMS.hpp Geodesic.hpp GeodesicLine.hpp \
	Gnomonic.hpp Math.hpp Utility.hpp
GeoidEval.o: GeoidEval.usage Config.h Constants.hpp DMS.hpp GeoCoords.hpp \
	Geoid.hpp Math.hpp UTMUPS.hpp Utility.hpp
MagneticField.o: MagneticField.usage CircularEngine.hpp Config.h Constants.hpp \
	DMS.hpp Geocentric.hpp MagneticModel.hpp Math.hpp SphericalEngine.hpp \
	SphericalHarmonic.hpp SphericalHarmonic1.hpp Utility.hpp
Planimeter.o: Planimeter.usage Accumulator.hpp Config.h Constants.hpp DMS.hpp \
	GeoCoords.hpp Geodesic.hpp Math.hpp PolygonArea.hpp UTMUPS.hpp \
	Utility.hpp
TransverseMercatorProj.o: TransverseMercatorProj.usage Config.h Constants.hpp \
	DMS.hpp EllipticFunction.hpp Math.hpp TransverseMercator.hpp \
	TransverseMercatorExact.hpp Utility.hpp

%: %.sh
	sed -e "s%@DEFAULTDIR@%$(GEOGRAPHICLIB_DATA)%" $< > $@
	chmod +x $@

INSTALL = install -b
DEST = $(PREFIX)/bin
SDEST = $(PREFIX)/sbin

install: $(PROGRAMS) $(SCRIPTS)
	test -f $(DEST) || mkdir -p $(DEST)
	$(INSTALL) $(PROGRAMS) $(DEST)
	test -f $(SDEST) || mkdir -p $(SDEST)
	$(INSTALL) $(SCRIPTS) $(SDEST)
