#
# Makefile.am
#
# Copyright (C) 2011, Charles Karney <charles@karney.com>
# $Id$

MATLAB_INTERFACE = \
geodesicdirect.cpp \
geodesicinverse.cpp \
geodesicline.cpp \
geoidheight.cpp \
mgrsforward.cpp \
mgrsreverse.cpp \
utmupsforward.cpp \
utmupsreverse.cpp \
geocentricforward.cpp \
geocentricreverse.cpp \
localcartesianforward.cpp \
localcartesianreverse.cpp \
polygonarea.cpp

MATLAB_DOC = \
geodesicdirect.m \
geodesicinverse.m \
geodesicline.m \
geoidheight.m \
mgrsforward.m \
mgrsreverse.m \
utmupsforward.m \
utmupsreverse.m \
geocentricforward.m \
geocentricreverse.m \
localcartesianforward.m \
localcartesianreverse.m \
polygonarea.m

MATLAB_COMPILESCRIPT = geographiclibinterface.m

MATLAB_ALL = $(MATLAB_INTERFACE) $(MATLAB_DOC) $(MATLAB_COMPILESCRIPT)

matlabdir=$(libexecdir)/GeographicLib/matlab

install:
	$(INSTALL) -d $(DESTDIR)$(matlabdir)
	$(INSTALL) -m 644 $(MATLAB_ALL) $(DESTDIR)$(matlabdir)

clean-local:
	rm -rf *.mex* *.oct

EXTRA_DIST = Makefile.mk CMakeLists.txt $(MATLAB_ALL)
