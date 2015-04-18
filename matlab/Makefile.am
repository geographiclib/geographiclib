#
# Makefile.am
#
# Copyright (C) 2011, Charles Karney <charles@karney.com>

MATLAB_FILES = \
$(srcdir)/geoddoc.m \
$(srcdir)/geodreckon.m \
$(srcdir)/geoddistance.m \
$(srcdir)/geodarea.m \
$(srcdir)/defaultellipsoid.m \
$(srcdir)/ecc2flat.m \
$(srcdir)/flat2ecc.m \
$(srcdir)/geodproj.m \
$(srcdir)/eqdazim_fwd.m \
$(srcdir)/eqdazim_inv.m \
$(srcdir)/cassini_fwd.m \
$(srcdir)/cassini_inv.m \
$(srcdir)/tranmerc_fwd.m \
$(srcdir)/tranmerc_inv.m \
$(srcdir)/gnomonic_fwd.m \
$(srcdir)/gnomonic_inv.m \
$(srcdir)/gedoc.m \
$(srcdir)/gereckon.m \
$(srcdir)/gedistance.m

MATLAB_LEGACY = \
$(srcdir)/geodesicdirect_a.m \
$(srcdir)/geodesicinverse_a.m \
$(srcdir)/geodesicline_a.m \
$(srcdir)/geoidheight_a.m \
$(srcdir)/mgrsforward_a.m \
$(srcdir)/mgrsreverse_a.m \
$(srcdir)/utmupsforward_a.m \
$(srcdir)/utmupsreverse_a.m \
$(srcdir)/geocentricforward_a.m \
$(srcdir)/geocentricreverse_a.m \
$(srcdir)/localcartesianforward_a.m \
$(srcdir)/localcartesianreverse_a.m \
$(srcdir)/polygonarea_a.m

MATLAB_PRIVATE = \
$(srcdir)/private/A1m1f.m \
$(srcdir)/private/A2m1f.m \
$(srcdir)/private/A3coeff.m \
$(srcdir)/private/A3f.m \
$(srcdir)/private/AngDiff.m \
$(srcdir)/private/AngNormalize2.m \
$(srcdir)/private/AngNormalize.m \
$(srcdir)/private/AngRound.m \
$(srcdir)/private/eatanhe.m \
$(srcdir)/private/C1f.m \
$(srcdir)/private/C1pf.m \
$(srcdir)/private/C2f.m \
$(srcdir)/private/C3coeff.m \
$(srcdir)/private/C3f.m \
$(srcdir)/private/C4coeff.m \
$(srcdir)/private/G4coeff.m \
$(srcdir)/private/C4f.m \
$(srcdir)/private/cbrtx.m \
$(srcdir)/private/cvmgt.m \
$(srcdir)/private/normalize.m \
$(srcdir)/private/SinCosSeries.m \
$(srcdir)/private/sumx.m \
$(srcdir)/private/swap.m

MATLAB_ALL = $(MATLAB_LEGACY) $(MATLAB_FILES)

matlabdir=$(libexecdir)/GeographicLib/matlab

install:
	$(INSTALL) -d $(DESTDIR)$(matlabdir)/private
	$(INSTALL) -m 644 $(MATLAB_ALL) $(DESTDIR)$(matlabdir)
	$(INSTALL) -m 644 $(MATLAB_PRIVATE) $(DESTDIR)$(matlabdir)/private

clean-local:
	rm -rf *.mex* *.oct

EXTRA_DIST = Makefile.mk CMakeLists.txt $(MATLAB_ALL) $(MATLAB_PRIVATE)
