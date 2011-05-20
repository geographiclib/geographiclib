# $Id$

LIBSTEM = Geographic
LIBRARY = lib$(LIBSTEM).a

all: $(LIBRARY)

INCLUDEPATH = ../include

PREFIX = /usr/local
GEOID_DEFAULT_PATH = $(PREFIX)/share/GeographicLib/geoids

MODULES = DMS EllipticFunction GeoCoords MGRS PolarStereographic \
	TransverseMercator TransverseMercatorExact UTMUPS Geocentric \
	LocalCartesian Geodesic GeodesicLine \
	AzimuthalEquidistant CassiniSoldner \
	Geoid LambertConformalConic Gnomonic OSGB AlbersEqualArea

HEADERS = Constants.hpp Config.h $(addsuffix .hpp,$(MODULES))
SOURCES = $(addsuffix .cpp,$(MODULES))
OBJECTS = $(addsuffix .o,$(MODULES))

CC = g++ -g
CXXFLAGS = -g -Wall -O3 -funroll-loops -finline-functions -fomit-frame-pointer

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

list:
	@echo $(SOURCES)
clean:
	rm -f *.o $(LIBRARY)

TAGS: $(HEADERS) $(SOURCES)
	etags $^

DMS.o: DMS.hpp Constants.hpp Config.h
EllipticFunction.o: EllipticFunction.hpp Constants.hpp Config.h
GeoCoords.o: GeoCoords.hpp Constants.hpp Config.h DMS.hpp MGRS.hpp UTMUPS.hpp
MGRS.o: MGRS.hpp Constants.hpp Config.h UTMUPS.hpp
PolarStereographic.o: PolarStereographic.hpp Constants.hpp Config.h
TransverseMercator.o: TransverseMercator.hpp Constants.hpp Config.h
TransverseMercatorExact.o: TransverseMercatorExact.hpp Constants.hpp Config.h \
	EllipticFunction.hpp
UTMUPS.o: UTMUPS.hpp Constants.hpp Config.h MGRS.hpp PolarStereographic.hpp \
	TransverseMercator.hpp
Geocentric.o: Geocentric.hpp Constants.hpp Config.h
LocalCartesian.o: LocalCartesian.hpp Constants.hpp Config.h Geocentric.hpp
Geodesic.o: Geodesic.hpp Constants.hpp Config.h GeodesicLine.hpp
GeodesicLine.o: GeodesicLine.hpp Constants.hpp Config.h Geodesic.hpp
AzimuthalEquidistant.o: AzimuthalEquidistant.hpp Constants.hpp Config.h \
	Geodesic.hpp
CassiniSoldner.o: CassiniSoldner.hpp Constants.hpp Config.h Geodesic.hpp
Geoid.o: Geoid.hpp Constants.hpp Config.h
LambertConformalConic.o: LambertConformalConic.hpp Constants.hpp Config.h
Gnomonic.o: Gnomonic.hpp Constants.hpp Config.h Geodesic.hpp
OSGB.o: OSGB.hpp Constants.hpp Config.h TransverseMercator.hpp
AlbersEqualArea.o: AlbersEqualArea.hpp Constants.hpp Config.h

.PHONY: all install list clean
