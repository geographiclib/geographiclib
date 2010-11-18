# $Id$

FUNCTIONS = utmupsforward utmupsreverse mgrsforward mgrsreverse \
	geodesicdirect geodesicinverse geodesicpath \
	geoidheight

MATLABFILES = $(addsuffix .cpp,$(FUNCTIONS)) $(addsuffix .m,$(FUNCTIONS))

all:
	@:
install:
	@:
clean:
	@:
list:
	@echo geographiclibinterface.m $(MATLABFILES)

.PHONY: all install list clean
