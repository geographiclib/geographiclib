# $Id$

FUNCTIONS = utmupsforward utmupsreverse mgrsforward mgrsreverse \
	geodesicdirect geodesicinverse geodesicline \
	geoidheight geocentricforward geocentricreverse \
	localcartesianforward localcartesianreverse

MATLABFILES = $(addsuffix .cpp,$(FUNCTIONS)) $(addsuffix .m,$(FUNCTIONS))

DEST = $(PREFIX)/share/GeographicLib/matlab
INSTALL=install -b

all:
	@:
install:
	test -d $(DEST) || mkdir -p $(DEST)
	$(INSTALL) -m 644 $(MATLABFILES) $(DEST)/
clean:
	@:
list:
	@echo geographiclibinterface.m $(MATLABFILES)

.PHONY: all install list clean
