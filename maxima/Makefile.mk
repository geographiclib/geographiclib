# $Id: 46ea48c30b46597075c4729396125e8bc95ba167 $

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
