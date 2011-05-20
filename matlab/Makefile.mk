# $Id: Makefile.mk 6896 2010-11-18 22:43:19Z karney $

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
