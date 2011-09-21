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

VPATH = ../include/GeographicLib

clean:
	rm -f *.o

ProjTest: ProjTest.o
TMTest: TMTest.o
GeodTest: GeodTest.o
ConicTest: ConicTest.o
NaNTester: NaNTester.o
PASouth: PASouth.o
ProjTest.o: Constants.hpp Math.hpp LambertConformalConic.hpp \
	PolarStereographic.hpp TransverseMercator.hpp \
	TransverseMercatorExact.hpp
TMTest.o: Constants.hpp Math.hpp TransverseMercator.hpp \
	TransverseMercatorExact.hpp Geodesic.hpp
GeodTest.o: Constants.hpp Math.hpp Geodesic.hpp
ConicTest.o: Constants.hpp Math.hpp LambertConformalConic.hpp \
	AlbersEqualArea.hpp
NaNTester.o: Constants.hpp Math.hpp EllipticFunction.hpp \
	TransverseMercator.hpp TransverseMercatorExact.hpp \
	PolarStereographic.hpp
PASouth.o: Constants.hpp Math.hpp LambertConformalConic.hpp DMS.hpp

INSTALL = install -b
PREFIX = /usr/local
DEST = $(PREFIX)/bin

install: $(PROGRAMS)
	test -f $(DEST) || mkdir -p $(DEST)
	$(INSTALL) $^ $(DEST)
