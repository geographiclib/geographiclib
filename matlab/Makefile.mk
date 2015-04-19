MATLAB_FILES = geoddoc.m geodreckon.m geoddistance.m geodarea.m \
	defaultellipsoid.m ecc2flat.m flat2ecc.m \
	geodproj.m eqdazim_fwd.m eqdazim_inv.m cassini_fwd.m cassini_inv.m \
	tranmerc_fwd.m tranmerc_inv.m gnomonic_fwd.m gnomonic_inv.m \
	gedoc.m gereckon.m gedistance.m

MATLAB_LEGACY = utmupsforward_a.m utmupsreverse_a.m mgrsforward_a.m \
	mgrsreverse_a.m geodesicdirect_a.m geodesicinverse_a.m \
	geodesicline_a.m geoidheight_a.m geocentricforward_a.m \
	geocentricreverse_a.m localcartesianforward_a.m \
	localcartesianreverse_a.m polygonarea_a.m

MATLAB_FILES = $(wildcard geographiclib/private/*.m)
MATLAB_PRIVATE = $(wildcard geographiclib/private/*.m)
MATLAB_LEGACY = $(wildcard geographiclib_legacy/*.m)

DEST = $(PREFIX)/libexec/matlab
INSTALL = install -b

all:
	@:

install:
	test -d $(DEST)/geographiclib/private || \
		mkdir -p $(DEST)/geographiclib/private
	test -d $(DEST)/geographiclib_legacy || \
		mkdir -p $(DEST)/geographiclib_legacy
	$(INSTALL) -m 644 $(MATLAB_FILES) $(DEST)/geographiclib
	$(INSTALL) -m 644 $(MATLAB_PRIVATE) $(DEST)/geographiclib/private/
	$(INSTALL) -m 644 $(MATLAB_LEGACY) $(DEST)/geographiclib_legacy

clean:
	rm -f *.mex* *.oct

.PHONY: all install clean
