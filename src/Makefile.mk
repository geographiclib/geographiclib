# $Id$

LIBSTEM = Geographic
LIBRARY = lib$(LIBSTEM).a

all: $(LIBRARY)

INCLUDEPATH = ../include

PREFIX = /usr/local
GEOID_DEFAULT_PATH = $(PREFIX)/share/GeographicLib/geoids

MODULES = AlbersEqualArea AzimuthalEquidistant CassiniSoldner \
	CircularEngine DMS EllipticFunction GeoCoords Geocentric \
	Geodesic GeodesicLine Geoid Gnomonic LambertConformalConic \
	LocalCartesian MGRS MagneticCircle MagneticModel OSGB \
	PolarStereographic PolygonArea SphericalEngine \
	TransverseMercator TransverseMercatorExact \
	UTMUPS
EXTRAHEADERS = Accumulator Constants Math SphericalHarmonic \
	SphericalHarmonic1 SphericalHarmonic2 Utility

HEADERS = Config.h $(addsuffix .hpp,$(EXTRAHEADERS) $(MODULES))
SOURCES = $(addsuffix .cpp,$(MODULES))
OBJECTS = $(addsuffix .o,$(MODULES))

CC = g++ -g
CXXFLAGS = -g -Wall -Wextra -O3

CPPFLAGS = -I$(INCLUDEPATH) $(DEFINES) \
	-DGEOID_DEFAULT_PATH=\"$(GEOID_DEFAULT_PATH)\"
LDFLAGS = $(LIBRARY)

$(LIBRARY): $(OBJECTS)
	$(AR) r $@ $?

VPATH = ../include/GeographicLib

INSTALL = install -b

install: $(LIBRARY)
	test -f $(PREFIX)/lib || mkdir -p $(PREFIX)/lib
	$(INSTALL) -m 644 $^ $(PREFIX)/lib

clean:
	rm -f *.o $(LIBRARY)

TAGS: $(HEADERS) $(SOURCES)
	etags $^

DMS.o: DMS.hpp Constants.hpp Math.hpp Utility.hpp Config.h
EllipticFunction.o: EllipticFunction.hpp Constants.hpp Math.hpp Config.h
GeoCoords.o: GeoCoords.hpp Constants.hpp Math.hpp Config.h DMS.hpp MGRS.hpp \
	UTMUPS.hpp
MGRS.o: MGRS.hpp Constants.hpp Math.hpp Utility.hpp Config.h UTMUPS.hpp
PolarStereographic.o: PolarStereographic.hpp Constants.hpp Math.hpp Config.h
TransverseMercator.o: TransverseMercator.hpp Constants.hpp Math.hpp Config.h
TransverseMercatorExact.o: TransverseMercatorExact.hpp Constants.hpp Math.hpp \
	Config.h EllipticFunction.hpp
UTMUPS.o: UTMUPS.hpp Constants.hpp Math.hpp Utility.hpp Config.h MGRS.hpp \
	PolarStereographic.hpp TransverseMercator.hpp
Geocentric.o: Geocentric.hpp Constants.hpp Math.hpp Config.h
LocalCartesian.o: LocalCartesian.hpp Constants.hpp Math.hpp Config.h \
	Geocentric.hpp
Geodesic.o: Geodesic.hpp Constants.hpp Math.hpp Config.h GeodesicLine.hpp
GeodesicLine.o: GeodesicLine.hpp Constants.hpp Math.hpp Config.h Geodesic.hpp
PolygonArea.o: PolygonArea.hpp Constants.hpp Math.hpp Config.h Geodesic.hpp \
	Accumulator.hpp
AzimuthalEquidistant.o: AzimuthalEquidistant.hpp Constants.hpp Math.hpp \
	Config.h Geodesic.hpp
CassiniSoldner.o: CassiniSoldner.hpp Constants.hpp Math.hpp Config.h \
	Geodesic.hpp
Geoid.o: Geoid.hpp Constants.hpp Math.hpp Config.h
LambertConformalConic.o: LambertConformalConic.hpp Constants.hpp Math.hpp \
	Config.h
Gnomonic.o: Gnomonic.hpp Constants.hpp Math.hpp Config.h Geodesic.hpp
OSGB.o: OSGB.hpp Constants.hpp Math.hpp Utility.hpp Config.h \
	TransverseMercator.hpp
AlbersEqualArea.o: AlbersEqualArea.hpp Constants.hpp Math.hpp Config.h
SphericalEngine.o: SphericalEngine.hpp Constants.hpp Math.hpp Config.h \
	CircularEngine.hpp
CircularEngine.o: CircularEngine.hpp Constants.hpp Math.hpp Config.h \
	SphericalEngine.hpp
MagneticModel.o: MagneticModel.hpp Constants.hpp Math.hpp Utility.hpp \
	Config.h SphericalEngine.hpp Geocentric.hpp MagneticCircle.hpp \
	SphericalHarmonic.hpp SphericalHarmonic1.hpp
MagneticCircle.o: MagneticCircle.hpp Constants.hpp Math.hpp Config.h \
	CircularEngine.hpp Geocentric.hpp

.PHONY: all install clean
