# $Id$

PROGRAMS = ProjTest TMTest

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

ProjTest: ProjTest.o
TMTest: TMTest.o

ProjTest.o: Constants.hpp LambertConformalConic.hpp PolarStereographic.hpp \
	TransverseMercator.hpp TransverseMercatorExact.hpp
TMTest.o: Constants.hpp TransverseMercator.hpp TransverseMercatorExact.hpp Geodesic.hpp

INSTALL = install -b
PREFIX = /usr/local
DEST = $(PREFIX)/bin

list:
	@echo $(addsuffix .cpp,$(PROGRAMS))

install: $(PROGRAMS)
	test -f $(DEST) || mkdir -p $(DEST)
	$(INSTALL) $^ $(DEST)

