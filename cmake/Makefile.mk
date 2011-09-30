# $Id: 9c58123e002564d0252e543c70539b5a79df36a6 $

DEST = $(PREFIX)/share/cmake/GeographicLib

INSTALL=install -b

all:
	@:
install:
	test -d $(DEST) || mkdir -p $(DEST)
	$(INSTALL) -m 644 FindGeographicLib.cmake $(DEST)
clean:
	@:

.PHONY: all install clean
