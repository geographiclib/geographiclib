# $Id$

FUNCTIONS = utmupsforward utmupsreverse mgrsforward mgrsreverse geoidheight

MATLABFILES = $(addsuffix .cpp,$(FUNCTIONS)) $(addsuffix .m,$(FUNCTIONS))

all:
	@:
install:
	@:
clean:
	@:
list:
	@echo compilematlabfuns.m $(MATLABFILES)

.PHONY: all install list clean
