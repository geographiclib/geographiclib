# $Id$
TARGET = GeoConvert TransverseMercatorTest ECEFConvert
all: $(TARGET)

CC = g++
CPPFLAGS = -I..
CXXFLAGS = -g -Wall -O0 -funroll-loops -finline-functions -fomit-frame-pointer

HEADERS = Constants.hpp DMS.hpp EllipticFunction.hpp GeoCoords.hpp MGRS.hpp \
	PolarStereographic.hpp TransverseMercator.hpp \
	TransverseMercatorExact.hpp UTMUPS.hpp ECEF.hpp LocalCartesian.hpp

SOURCES = $(patsubst %.hpp,%.cpp,$(HEADERS)) \
	GeoConvert.cpp TransverseMercatorTest.cpp ECEFConvert.cpp

GeoConvert: GeoConvert.o GeoCoords.o MGRS.o UTMUPS.o DMS.o Constants.o \
	TransverseMercator.o PolarStereographic.o
TransverseMercatorTest: TransverseMercatorTest.o TransverseMercatorExact.o \
	Constants.o EllipticFunction.o TransverseMercator.o
ECEFConvert: ECEFConvert.o ECEF.o LocalCartesian.o Constants.o

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
ECEF.o: ECEF.hpp Constants.hpp
LocalCartesian.o: LocalCartesian.hpp ECEF.hpp Constants.hpp
ECEFConvert.o: ECEF.hpp LocalCartesian.hpp Constants.hpp
GeoConvert.o: GeoCoords.hpp UTMUPS.hpp
TransverseMercatorTest.o: TransverseMercatorExact.hpp EllipticFunction.hpp \
	TransverseMercator.hpp Constants.hpp

FIGURES = gauss-krueger-graticule thompson-tm-graticule \
	gauss-krueger-convergence-scale gauss-schreiber-graticule-a \
	gauss-krueger-graticule-a thompson-tm-graticule-a

MAXIMASOURCES = tm.mac ellint.mac tmseries.mac

doc: Doxyfile Geographic.doc \
	$(HEADERS) $(SOURCES) \
	$(addsuffix .pdf,$(FIGURES)) $(addsuffix .png,$(FIGURES))
	rm -rf html/*
	doxygen
	for f in $(FIGURES); do cp -p $$f.pdf $$f.png html/;done
	for f in $(MAXIMASOURCES); do cp -p $$f html/;done
	touch $@

clean:
	rm -f *.o
