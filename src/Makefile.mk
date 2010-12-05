# $Id$

LIBSTEM = Geographic
LIBRARY = lib$(LIBSTEM).a

all: $(LIBRARY)

INCLUDEPATH = ../include

MODULES = DMS EllipticFunction GeoCoords MGRS PolarStereographic \
	TransverseMercator TransverseMercatorExact UTMUPS Geocentric \
	LocalCartesian Geodesic GeodesicLine \
	AzimuthalEquidistant CassiniSoldner \
	Geoid LambertConformalConic Gnomonic OSGB AlbersEqualArea

HEADERS = Constants.hpp $(addsuffix .hpp,$(MODULES))
SOURCES = $(addsuffix .cpp,$(MODULES))
OBJECTS = $(addsuffix .o,$(MODULES))

CC = g++ -g
CXXFLAGS = -g -Wall -O3 -funroll-loops -finline-functions -fomit-frame-pointer

CPPFLAGS = -I$(INCLUDEPATH) $(DEFINES)
LDFLAGS = $(LIBRARY)

$(LIBRARY): $(OBJECTS)
	$(AR) r $@ $?

VPATH = ../include/GeographicLib

INSTALL = install -b
PREFIX = /usr/local

install: $(LIBRARY)
	test -f $(PREFIX)/lib || mkdir -p $(PREFIX)/lib
	$(INSTALL) -m 644 $^ $(PREFIX)/lib

list:
	@echo $(SOURCES)
clean:
	rm -f *.o $(LIBRARY)

TAGS: $(HEADERS) $(SOURCES)
	etags $^

DMS.o: DMS.hpp Constants.hpp
EllipticFunction.o: EllipticFunction.hpp Constants.hpp
GeoCoords.o: GeoCoords.hpp Constants.hpp DMS.hpp MGRS.hpp UTMUPS.hpp
MGRS.o: MGRS.hpp Constants.hpp UTMUPS.hpp
PolarStereographic.o: PolarStereographic.hpp Constants.hpp
TransverseMercator.o: TransverseMercator.hpp Constants.hpp
TransverseMercatorExact.o: TransverseMercatorExact.hpp Constants.hpp \
	EllipticFunction.hpp
UTMUPS.o: UTMUPS.hpp Constants.hpp MGRS.hpp PolarStereographic.hpp \
	TransverseMercator.hpp
Geocentric.o: Geocentric.hpp Constants.hpp
LocalCartesian.o: LocalCartesian.hpp Constants.hpp Geocentric.hpp
Geodesic.o: Geodesic.hpp Constants.hpp GeodesicLine.hpp
GeodesicLine.o: GeodesicLine.hpp Constants.hpp Geodesic.hpp
AzimuthalEquidistant.o: AzimuthalEquidistant.hpp Constants.hpp Geodesic.hpp
CassiniSoldner.o: CassiniSoldner.hpp Constants.hpp Geodesic.hpp
Geoid.o: Geoid.hpp Constants.hpp
LambertConformalConic.o: LambertConformalConic.hpp Constants.hpp
Gnomonic.o: Gnomonic.hpp Constants.hpp Geodesic.hpp
OSGB.o: OSGB.hpp Constants.hpp TransverseMercator.hpp
AlbersEqualArea.o: AlbersEqualArea.hpp Constants.hpp

.PHONY: all install list clean
