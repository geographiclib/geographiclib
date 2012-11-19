FUNCTIONS = utmupsforward utmupsreverse mgrsforward mgrsreverse \
	geodesicdirect geodesicinverse geodesicline \
	geoidheight geocentricforward geocentricreverse \
	localcartesianforward localcartesianreverse polygonarea

MATLAB_COMPILESCRIPT = geographiclibinterface.m

MATLAB_GEOD = geoddoc.m geodreckon.m geoddistance.m geodarea.m \
	defaultellipsoid.m
MATLAB_GEOD_PRIVATE = $(wildcard private/*.m)

MATLABFILES = $(addsuffix .cpp,$(FUNCTIONS)) $(addsuffix .m,$(FUNCTIONS)) \
	 $(MATLAB_COMPILESCRIPT) $(MATLAB_GEOD)

DEST = $(PREFIX)/libexec/GeographicLib/matlab
INSTALL = install -b

all:
	@:

install:
	test -d $(DEST)/private || mkdir -p $(DEST)/private
	$(INSTALL) -m 644 $(MATLABFILES) $(DEST)/
	$(INSTALL) -m 644 $(MATLAB_GEOD_PRIVATE) $(DEST)/private/
clean:
	rm -f *.mex* *.oct

.PHONY: all install clean
