# $Id: c26d5af8f3ff807ccf3f27ec7627cf21d1efa5c9 $

FUNCTIONS = utmupsforward utmupsreverse mgrsforward mgrsreverse \
	geodesicdirect geodesicinverse geodesicline \
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
