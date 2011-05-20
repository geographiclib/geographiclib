# $Id: Makefile.mk 6714 2009-10-17 11:50:59Z ckarney $

MAXIMA = tm ellint tmseries geod
MAXIMASOURCES = $(addsuffix .mac,$(MAXIMA))

all:
	@:
install:
	@:
clean:
	@:
list:
	@echo $(MAXIMASOURCES)

.PHONY: all install list clean
