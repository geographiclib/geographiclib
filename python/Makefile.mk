# $Id: c57d9935a4847a2f5564fc68da6b3275ba36221d $

MODULES = __init__ geomath constants accumulator geodesiccapability \
	geodesic geodesicline polygonarea
PACKAGE=geographiclib
PYTHON_FILES = $(patsubst %,$(PACKAGE)/%.py,$(MODULES))

DEST = $(PREFIX)/lib/python/site-packages/$(PACKAGE)
INSTALL = install -b

all:
	@:

install:
	test -d $(DEST) || mkdir -p $(DEST)
	$(INSTALL) -m 644 $(PYTHON_FILES) $(DEST)/

clean:
	rm -f *.pyc $(PACKAGE)/*.pyc

list:
	@echo Makefile.am $(PACKAGE)/CMakeLists.txt $(PYTHONFILES)

.PHONY: all install list clean
