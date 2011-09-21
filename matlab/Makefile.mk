# $Id$

FUNCTIONS = utmupsforward utmupsreverse mgrsforward mgrsreverse \
	geodesicdirect geodesicinverse geodesicline \
	geoidheight geocentricforward geocentricreverse \
	localcartesianforward localcartesianreverse polygonarea

MATLAB_COMPILESCRIPT = geographiclibinterface.m

MATLABFILES = $(addsuffix .cpp,$(FUNCTIONS)) $(addsuffix .m,$(FUNCTIONS)) \
	 $(MATLAB_COMPILESCRIPT)

DEST = $(PREFIX)/libexec/GeographicLib/matlab
INSTALL = install -b

all:
	@:

install:
	test -d $(DEST) || mkdir -p $(DEST)
	$(INSTALL) -m 644 $(MATLABFILES) $(DEST)/
clean:
	rm -f *.mex* *.oct

.PHONY: all install clean
