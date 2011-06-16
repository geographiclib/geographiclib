# $Id$

PROGRAMS = ProjTest TMTest GeodTest ConicTest NaNTester PASouth

all: $(PROGRAMS)

LIBSTEM = Geographic
LIBRARY = lib$(LIBSTEM).a

INCLUDEPATH = ../include
LIBPATH = ../src

# After installation, use these values of INCLUDEPATH and LIBPATH
# INCLUDEPATH = $(PREFIX)/include
# LIBPATH = $(PREFIX)/lib

CC = g++ -g
CXXFLAGS = -g -Wall -Wextra -O3

CPPFLAGS = -I$(INCLUDEPATH) $(DEFINES)
LDLIBS = -L$(LIBPATH) -l$(LIBSTEM)

$(PROGRAMS): $(LIBPATH)/$(LIBRARY)
	$(CC) -o $@ $@.o $(LDLIBS)

GeodTestL: $(LIBPATH)/$(LIBRARY) ../srcL/libGeographicL.a
	$(CC) -o $@ $@.o $(LDLIBS) -L$(LIBPATH)L -l$(LIBSTEM)L

VPATH = ../include/GeographicLib

clean:
	rm -f *.o

ProjTest: ProjTest.o
TMTest: TMTest.o
GeodTest: GeodTest.o
GeodTestL: GeodTestL.o
ConicTest: ConicTest.o
NaNTester: NaNTester.o
PASouth: PASouth.o
ProjTest.o: Constants.hpp LambertConformalConic.hpp PolarStereographic.hpp \
	TransverseMercator.hpp TransverseMercatorExact.hpp
TMTest.o: Constants.hpp TransverseMercator.hpp TransverseMercatorExact.hpp \
	Geodesic.hpp
GeodTest.o: Constants.hpp Geodesic.hpp
ConicTest.o: Constants.hpp LambertConformalConic.hpp AlbersEqualArea.hpp
NaNTester.o: Constants.hpp EllipticFunction.hpp TransverseMercator.hpp \
	TransverseMercatorExact.hpp PolarStereographic.hpp
PASouth.o: Constants.hpp LambertConformalConic.hpp DMS.hpp

GeodTestL.o: GeodTest.cpp Constants.hpp Geodesic.hpp \
	../include/GeographicLibL/Constants.hpp \
	../include/GeographicLibL/Geodesic.hpp
	$(CC) $(CXXFLAGS) -I$(INCLUDEPATH) -DUSE_LONG_DOUBLE_GEOGRAPHICLIB \
	-c -o $@ $<

INSTALL = install -b
PREFIX = /usr/local
DEST = $(PREFIX)/bin

list:
	@echo $(addsuffix .cpp,$(PROGRAMS))

install: $(PROGRAMS)
	test -f $(DEST) || mkdir -p $(DEST)
	$(INSTALL) $^ $(DEST)

