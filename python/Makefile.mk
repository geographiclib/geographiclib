# $Id: a8536304680ae3fee858ff4a5c78a081ea9a049b $

PYTHONFILES = GeographicLib.py

DEST = $(PREFIX)/share/GeographicLib/python
INSTALL = install -b

all:
	@:

install:
	test -d $(DEST) || mkdir -p $(DEST)
	$(INSTALL) -m 644 $(PYTHONFILES) $(DEST)/

clean:
	rm -f *.pyc

list:
	@echo $(PYTHONFILES)

.PHONY: all install list clean
