MATLAB_FILES = $(wildcard geographiclib/*.m)
MATLAB_PRIVATE = $(wildcard geographiclib/private/*.m)
MATLAB_LEGACY = $(wildcard geographiclib_legacy/*.m) \
	$(wildcard geographiclib_legacy/*.cpp)

DEST = $(PREFIX)/share/matlab
INSTALL = install -b

all:
	@:

install:
	test -d $(DEST)/geographiclib/private || \
		mkdir -p $(DEST)/geographiclib/private
	test -d $(DEST)/geographiclib_legacy || \
		mkdir -p $(DEST)/geographiclib_legacy
	$(INSTALL) -m 644 $(MATLAB_FILES) $(DEST)/geographiclib
	$(INSTALL) -m 644 $(MATLAB_PRIVATE) $(DEST)/geographiclib/private/
	$(INSTALL) -m 644 $(MATLAB_LEGACY) $(DEST)/geographiclib_legacy

clean:
	rm -f *.mex* *.oct

.PHONY: all install clean
