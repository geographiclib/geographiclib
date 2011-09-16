# $Id$

VERSION = 9.0.1

TEMPLATE = lib

INCLUDEPATH = ../include
INCLUDEDIR = $$INCLUDEPATH/GeographicLib

SOURCES += DMS.cpp
SOURCES += EllipticFunction.cpp
SOURCES += GeoCoords.cpp
SOURCES += MGRS.cpp
SOURCES += PolarStereographic.cpp
SOURCES += TransverseMercator.cpp
SOURCES += TransverseMercatorExact.cpp
SOURCES += UTMUPS.cpp
SOURCES += Geocentric.cpp
SOURCES += LocalCartesian.cpp
SOURCES += Geodesic.cpp
SOURCES += GeodesicLine.cpp
SOURCES += PolygonArea.cpp
SOURCES += AzimuthalEquidistant.cpp
SOURCES += CassiniSoldner.cpp
SOURCES += Geoid.cpp
SOURCES += LambertConformalConic.cpp
SOURCES += Gnomonic.cpp
SOURCES += OSGB.cpp
SOURCES += AlbersEqualArea.cpp

HEADERS += $$INCLUDEDIR/DMS.hpp
HEADERS += $$INCLUDEDIR/EllipticFunction.hpp
HEADERS += $$INCLUDEDIR/GeoCoords.hpp
HEADERS += $$INCLUDEDIR/MGRS.hpp
HEADERS += $$INCLUDEDIR/PolarStereographic.hpp
HEADERS += $$INCLUDEDIR/TransverseMercator.hpp
HEADERS += $$INCLUDEDIR/TransverseMercatorExact.hpp
HEADERS += $$INCLUDEDIR/UTMUPS.hpp
HEADERS += $$INCLUDEDIR/Geocentric.hpp
HEADERS += $$INCLUDEDIR/LocalCartesian.hpp
HEADERS += $$INCLUDEDIR/Geodesic.hpp
HEADERS += $$INCLUDEDIR/GeodesicLine.hpp
HEADERS += $$INCLUDEDIR/PolygonArea.hpp
HEADERS += $$INCLUDEDIR/AzimuthalEquidistant.hpp
HEADERS += $$INCLUDEDIR/CassiniSoldner.hpp
HEADERS += $$INCLUDEDIR/Geoid.hpp
HEADERS += $$INCLUDEDIR/LambertConformalConic.hpp
HEADERS += $$INCLUDEDIR/Gnomonic.hpp
HEADERS += $$INCLUDEDIR/OSGB.hpp
HEADERS += $$INCLUDEDIR/AlbersEqualArea.hpp

HEADERS += $$INCLUDEDIR/Config.h
HEADERS += $$INCLUDEDIR/Constants.hpp
HEADERS += $$INCLUDEDIR/Math.hpp
HEADERS += $$INCLUDEDIR/Accumulator.hpp
