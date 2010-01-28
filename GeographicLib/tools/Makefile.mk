# $Id$

PROGRAMS = GeoConvert TransverseMercatorTest CartConvert Geod EquidistantTest \
	GeoidEval

all: $(PROGRAMS)

LIBSTEM = Geographic
LIBRARY = lib$(LIBSTEM).a

INCLUDEPATH = ../include
LIBPATH = ../src

# After installation, use these values of INCLUDEPATH and LIBPATH
# INCLUDEPATH = $(PREFIX)/include
# LIBPATH = $(PREFIX)/lib

CC = g++ -g
CXXFLAGS = -g -Wall -O3 -funroll-loops -finline-functions -fomit-frame-pointer

CPPFLAGS = -I$(INCLUDEPATH) $(DEFINES)
LDLIBS = -L$(LIBPATH) -l$(LIBSTEM)


$(PROGRAMS): $(LIBPATH)/$(LIBRARY)
	$(CC) -o $@ $@.o $(LDLIBS)

VPATH = ../include/GeographicLib

clean:
	rm -f *.o

GeoConvert: GeoConvert.o
TransverseMercatorTest: TransverseMercatorTest.o
CartConvert: CartConvert.o
Geod: Geod.o
EquidistantTest: EquidistantTest.o
GeoidEval: GeoidEval.o

GeoConvert.o: Constants.hpp GeoCoords.hpp UTMUPS.hpp
TransverseMercatorTest.o: Constants.hpp EllipticFunction.hpp \
	TransverseMercator.hpp TransverseMercatorExact.hpp
CartConvert.o: Constants.hpp Geocentric.hpp LocalCartesian.hpp
Geod.o: Constants.hpp DMS.hpp Geodesic.hpp
EquidistantTest.o: AzimuthalEquidistant.hpp CassiniSoldner.hpp Constants.hpp \
	Geodesic.hpp
GeoidEval.o: Constants.hpp DMS.hpp Geoid.hpp

INSTALL = install -b
PREFIX = /usr/local
DEST = $(PREFIX)/bin

list:
	@echo $(addsuffix .cpp,$(PROGRAMS))

install: $(PROGRAMS)
	test -f $(DEST) || mkdir -p $(DEST)
	$(INSTALL) $^ $(DEST)

