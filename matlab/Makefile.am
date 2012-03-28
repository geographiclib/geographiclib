#
# Makefile.am
#
# Copyright (C) 2011, Charles Karney <charles@karney.com>
# $Id$

MATLAB_INTERFACE = \
$(srcdir)/geodesicdirect.cpp \
$(srcdir)/geodesicinverse.cpp \
$(srcdir)/geodesicline.cpp \
$(srcdir)/geoidheight.cpp \
$(srcdir)/mgrsforward.cpp \
$(srcdir)/mgrsreverse.cpp \
$(srcdir)/utmupsforward.cpp \
$(srcdir)/utmupsreverse.cpp \
$(srcdir)/geocentricforward.cpp \
$(srcdir)/geocentricreverse.cpp \
$(srcdir)/localcartesianforward.cpp \
$(srcdir)/localcartesianreverse.cpp \
$(srcdir)/polygonarea.cpp

MATLAB_DOC = \
$(srcdir)/geodesicdirect.m \
$(srcdir)/geodesicinverse.m \
$(srcdir)/geodesicline.m \
$(srcdir)/geoidheight.m \
$(srcdir)/mgrsforward.m \
$(srcdir)/mgrsreverse.m \
$(srcdir)/utmupsforward.m \
$(srcdir)/utmupsreverse.m \
$(srcdir)/geocentricforward.m \
$(srcdir)/geocentricreverse.m \
$(srcdir)/localcartesianforward.m \
$(srcdir)/localcartesianreverse.m \
$(srcdir)/polygonarea.m

MATLAB_COMPILESCRIPT = $(srcdir)/geographiclibinterface.m

MATLAB_ALL = $(MATLAB_INTERFACE) $(MATLAB_DOC) $(MATLAB_COMPILESCRIPT)

matlabdir=$(libexecdir)/GeographicLib/matlab

install:
	$(INSTALL) -d $(DESTDIR)$(matlabdir)
	$(INSTALL) -m 644 $(MATLAB_ALL) $(DESTDIR)$(matlabdir)

clean-local:
	rm -rf *.mex* *.oct

EXTRA_DIST = Makefile.mk CMakeLists.txt $(MATLAB_ALL)
