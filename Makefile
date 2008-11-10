# $Id$
TARGET = GeoConvert TransverseMercatorTest
all: $(TARGET)

CC = g++
CPPFLAGS = -I..
CXXFLAGS = -g -Wall -O3 -funroll-loops -finline-functions -fomit-frame-pointer

GeoConvert: GeoConvert.o GeoCoords.o MGRS.o UTMUPS.o DMS.o Constants.o \
	TransverseMercator.o PolarStereographic.o
TransverseMercatorTest: TransverseMercatorTest.o TransverseMercatorExact.o \
	Constants.o EllipticFunction.o

Constants.o: Constants.hpp
DMS.o: DMS.hpp
EllipticFunction.o: EllipticFunction.hpp Constants.hpp
GeoConvert.o: GeoCoords.hpp UTMUPS.hpp
GeoCoords.o: GeoCoords.hpp UTMUPS.hpp MGRS.hpp DMS.hpp
MGRS.o: MGRS.hpp UTMUPS.hpp
PolarStereographic.o: PolarStereographic.hpp Constants.hpp
TransverseMercator.o: TransverseMercator.hpp Constants.hpp
TransverseMercatorExact.o: TransverseMercatorExact.hpp EllipticFunction.hpp Constants.hpp
TransverseMercatorTest.o: TransverseMercatorExact.hpp EllipticFunction.hpp Constants.hpp
UTMUPS.o: UTMUPS.hpp MGRS.hpp PolarStereographic.hpp TransverseMercator.hpp

doc: Doxyfile Geographic.doc \
	GeoConvert.cpp TransverseMercatorTest.cpp \
	Constants.cpp DMS.cpp EllipticFunction.cpp GeoCoords.cpp MGRS.cpp \
	PolarStereographic.cpp TransverseMercator.cpp TransverseMercatorExact.cpp UTMUPS.cpp \
	Constants.hpp DMS.hpp EllipticFunction.hpp GeoCoords.hpp MGRS.hpp \
	PolarStereographic.hpp TransverseMercator.hpp TransverseMercatorExact.hpp UTMUPS.hpp
	rm -rf html/*
	doxygen
	touch $@

clean:
	rm -f *.o

# LDLIBS = -lproj

