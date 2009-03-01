# $Id$

LIBSTEM = Geographic
LIBRARY = lib$(LIBSTEM).a
PROGRAMS = GeoConvert TransverseMercatorTest CartConvert Geod

all: $(PROGRAMS) $(LIBRARY)

CC = g++ -g
CXXFLAGS = -g -Wall -O3 -funroll-loops -finline-functions -fomit-frame-pointer

CPPFLAGS = -I..
LDFLAGS = $(LIBRARY)

PREFIX = /usr/local
# After installation, use these values of CPPFLAGS and LDFLAGS

# CPPFLAGS = -I$(PREFIX)/include
# LDFLAGS = -L$(PREFIX)/lib -l$(LIBSTEM)

MODULES = DMS EllipticFunction GeoCoords MGRS PolarStereographic \
	TransverseMercator TransverseMercatorExact UTMUPS Geocentric \
	LocalCartesian Geodesic

HEADERS = Constants.hpp $(patsubst %,%.hpp,$(MODULES))
SOURCES = $(patsubst %,%.cpp,$(MODULES))
OBJECTS = $(patsubst %,%.o,$(MODULES))
ALLSOURCES = $(SOURCES) $(patsubst %,%.cpp,$(PROGRAMS))

$(LIBRARY): $(OBJECTS)
	$(AR) r $@ $?

$(PROGRAMS): $(LIBRARY)
	$(CC) -o $@ $@.o $(LDFLAGS)

GeoConvert: GeoConvert.o
TransverseMercatorTest: TransverseMercatorTest.o
CartConvert: CartConvert.o
Geod: Geod.o

Constants.o: Constants.hpp
DMS.o: DMS.hpp
EllipticFunction.o: EllipticFunction.hpp Constants.hpp
GeoCoords.o: GeoCoords.hpp UTMUPS.hpp MGRS.hpp DMS.hpp
MGRS.o: MGRS.hpp UTMUPS.hpp
PolarStereographic.o: PolarStereographic.hpp Constants.hpp
TransverseMercator.o: TransverseMercator.hpp Constants.hpp
TransverseMercatorExact.o: TransverseMercatorExact.hpp EllipticFunction.hpp \
	Constants.hpp
UTMUPS.o: UTMUPS.hpp MGRS.hpp PolarStereographic.hpp TransverseMercator.hpp
Geocentric.o: Geocentric.hpp Constants.hpp
LocalCartesian.o: LocalCartesian.hpp Geocentric.hpp Constants.hpp
Geodesic.o: Geodesic.hpp Constants.hpp
GeoConvert.o: GeoCoords.hpp UTMUPS.hpp
TransverseMercatorTest.o: EllipticFunction.hpp TransverseMercatorExact.hpp \
	TransverseMercator.hpp
CartConvert.o: Geocentric.hpp LocalCartesian.hpp
Geod.o: Geodesic.hpp DMS.hpp

FIGURES = gauss-krueger-graticule thompson-tm-graticule \
	gauss-krueger-convergence-scale gauss-schreiber-graticule-a \
	gauss-krueger-graticule-a thompson-tm-graticule-a

MAXIMASOURCES = tm.mac ellint.mac tmseries.mac

install: install-lib install-headers install-progs

install-lib: $(LIBRARY)
	test -f $(PREFIX)/lib || mkdir -p $(PREFIX)/lib
	install -m 644 $(LIBRARY) $(PREFIX)/lib

install-headers: $(HEADERS)
	test -f $(PREFIX)/include/GeographicLib || mkdir -p $(PREFIX)/include/GeographicLib
	install -m 644 $(HEADERS) $(PREFIX)/include/GeographicLib

install-progs: $(PROGRAMS)
	test -f $(PREFIX)/bin || mkdir -p $(PREFIX)/bin
	install $(PROGRAMS) $(PREFIX)/bin

doc: Doxyfile Geographic.doc \
	$(HEADERS) $(ALLSOURCES) \
	$(addsuffix .pdf,$(FIGURES)) $(addsuffix .png,$(FIGURES))
	rm -rf html/*
	doxygen
	for f in $(FIGURES); do cp -p $$f.pdf $$f.png html/;done
	for f in $(MAXIMASOURCES); do cp -p $$f html/;done
	touch $@

clean:
	rm -f *.o $(LIBRARY)
