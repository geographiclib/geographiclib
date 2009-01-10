# $Id$
TARGET = GeoConvert TransverseMercatorTest
all: $(TARGET)

CC = g++
CPPFLAGS = -I..
CXXFLAGS = -g -Wall -O3 -funroll-loops -finline-functions -fomit-frame-pointer


HEADERS = Constants.hpp DMS.hpp EllipticFunction.hpp GeoCoords.hpp MGRS.hpp \
	PolarStereographic.hpp TransverseMercator.hpp \
	TransverseMercatorExact.hpp UTMUPS.hpp

SOURCES = Constants.cpp DMS.cpp EllipticFunction.cpp GeoCoords.cpp MGRS.cpp \
	PolarStereographic.cpp TransverseMercator.cpp \
	TransverseMercatorExact.cpp UTMUPS.cpp \
	GeoConvert.cpp TransverseMercatorTest.cpp

GeoConvert: GeoConvert.o GeoCoords.o MGRS.o UTMUPS.o DMS.o Constants.o \
	TransverseMercator.o PolarStereographic.o
TransverseMercatorTest: TransverseMercatorTest.o TransverseMercatorExact.o \
	Constants.o EllipticFunction.o TransverseMercator.o

Constants.o: Constants.hpp
DMS.o: DMS.hpp
EllipticFunction.o: EllipticFunction.hpp Constants.hpp
GeoConvert.o: GeoCoords.hpp UTMUPS.hpp
GeoCoords.o: GeoCoords.hpp UTMUPS.hpp MGRS.hpp DMS.hpp
MGRS.o: MGRS.hpp UTMUPS.hpp
PolarStereographic.o: PolarStereographic.hpp Constants.hpp
TransverseMercator.o: TransverseMercator.hpp Constants.hpp
TransverseMercatorExact.o: TransverseMercatorExact.hpp EllipticFunction.hpp \
	Constants.hpp
TransverseMercatorTest.o: TransverseMercatorExact.hpp EllipticFunction.hpp \
	Constants.hpp TransverseMercator.hpp
UTMUPS.o: UTMUPS.hpp MGRS.hpp PolarStereographic.hpp TransverseMercator.hpp

FIGURES = gauss-krueger-graticule thompson-tm-graticule \
	gauss-krueger-convergence-scale gauss-laborde-graticule-a \
	gauss-krueger-graticule-a thompson-tm-graticule-a

doc: Doxyfile Geographic.doc \
	$(HEADERS) $(SOURCES) \
	$(addsuffix .pdf,$(FIGURES)) $(addsuffix .png,$(FIGURES))
	rm -rf html/*
	doxygen
	for f in $(FIGURES); do cp $$f.pdf $$f.png html/;done
	touch $@

clean:
	rm -f *.o
