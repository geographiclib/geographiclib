# $Id: Makefile.mk 6859 2010-09-06 14:45:33Z karney $

LIBSTEM = Geographic
LIBRARY = lib$(LIBSTEM).a

all: $(LIBRARY)

INCLUDEPATH = ../include

MODULES = DMS EllipticFunction GeoCoords MGRS PolarStereographic \
	TransverseMercator TransverseMercatorExact UTMUPS Geocentric \
	LocalCartesian Geodesic GeodesicLine \
	AzimuthalEquidistant CassiniSoldner \
	Geoid LambertConformalConic Gnomonic

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

DMS.o: Constants.hpp DMS.hpp
EllipticFunction.o: Constants.hpp EllipticFunction.hpp
GeoCoords.o: Constants.hpp DMS.hpp GeoCoords.hpp MGRS.hpp UTMUPS.hpp
MGRS.o: Constants.hpp MGRS.hpp UTMUPS.hpp
PolarStereographic.o: Constants.hpp PolarStereographic.hpp
TransverseMercator.o: Constants.hpp TransverseMercator.hpp
TransverseMercatorExact.o: Constants.hpp EllipticFunction.hpp \
	TransverseMercatorExact.hpp
UTMUPS.o: Constants.hpp MGRS.hpp PolarStereographic.hpp \
	TransverseMercator.hpp UTMUPS.hpp
Geocentric.o: Constants.hpp Geocentric.hpp
LocalCartesian.o: Constants.hpp Geocentric.hpp LocalCartesian.hpp
Geodesic.o: Constants.hpp Geodesic.hpp GeodesicLine.hpp
GeodesicLine.o: Constants.hpp Geodesic.hpp GeodesicLine.hpp
AzimuthalEquidistant.o: AzimuthalEquidistant.hpp Constants.hpp Geodesic.hpp
CassiniSoldner.o: CassiniSoldner.hpp Constants.hpp Geodesic.hpp
Geoid.o: Constants.hpp Geoid.hpp
LambertConformalConic.o: Constants.hpp LambertConformalConic.hpp
Gnomonic.o: Gnomonic.hpp Constants.hpp Geodesic.hpp

.PHONY: all install list clean
