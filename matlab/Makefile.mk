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

MATLAB_PRIVATE = $(wildcard private/*.m)

MATLAB_ALL = $(MATLAB_FILES) $(MATLAB_LEGACY)

DEST = $(PREFIX)/libexec/GeographicLib/matlab
INSTALL = install -b

all:
	@:

install:
	test -d $(DEST)/private || mkdir -p $(DEST)/private
	$(INSTALL) -m 644 $(MATLAB_ALL) $(DEST)/
	$(INSTALL) -m 644 $(MATLAB_PRIVATE) $(DEST)/private/
clean:
	rm -f *.mex* *.oct

.PHONY: all install clean
