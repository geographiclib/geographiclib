#
# Makefile.am
#
# Copyright (C) 2009, Francesco P. Lovergine <frankie@debian.org>
# $Id$

nobase_include_HEADERS = GeographicLib/AzimuthalEquidistant.hpp \
				GeographicLib/CassiniSoldner.hpp \
				GeographicLib/Constants.hpp \
				GeographicLib/DMS.hpp \
				GeographicLib/EllipticFunction.hpp \
				GeographicLib/Geocentric.hpp \
				GeographicLib/GeoCoords.hpp \
				GeographicLib/Geodesic.hpp \
				GeographicLib/GeodesicLine.hpp \
				GeographicLib/Geoid.hpp \
				GeographicLib/Gnomonic.hpp \
				GeographicLib/LambertConformalConic.hpp \
				GeographicLib/LocalCartesian.hpp \
				GeographicLib/MGRS.hpp \
				GeographicLib/OSGB.hpp \
				GeographicLib/PolarStereographic.hpp \
				GeographicLib/TransverseMercatorExact.hpp \
				GeographicLib/TransverseMercator.hpp \
				GeographicLib/UTMUPS.hpp

geoids_dir=$(datadir)/GeographicLib/geoids

DEFS=-DGEOID_DEFAULT_PATH=\"$(geoids_dir)\"  @DEFS@

EXTRA_DIST = Makefile.mk
