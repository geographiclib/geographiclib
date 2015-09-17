JS_VERSION := $(shell grep "version_string = " Math.js | cut -f2 -d'"')
SCRIPTS = Math.js Geodesic.js GeodesicLine.js PolygonArea.js DMS.js
TUTORIALS = tutorials/geodesics.md tutorials/interface.md \
	tutorials/examples.md tutorials/tutorials.json
ALLSCRIPTS = $(DEST)/geographiclib.js $(DEST)/geographiclib.min.js \
	$(patsubst %,$(DEST)/src/%,$(SCRIPTS))
DEST=node_modules/geographiclib
NODE_VERSION := $(shell grep "version.: " $(DEST)/package.json | cut -f4 -d'"')
NODEFILES = \
	$(DEST)/package.json \
	$(DEST)/README.md \
	$(DEST)/LICENSE.txt \
	$(ALLSCRIPTS)

all: doc script samples test pack
doc: html/index.html
script: $(ALLSCRIPTS)
pack: geographiclib-$(NODE_VERSION).tgz

# requires that jsdoc be installed
html/index.html: GeographicLib.md $(SCRIPTS) $(TUTORIALS)
	jsdoc -d html -u tutorials -R GeographicLib.md $(SCRIPTS)

publish: doc
	rsync --delete -a -e ssh html/ karney,geographiclib@web.sourceforge.net:htdocs/$(JS_VERSION)/js/

geographiclib-$(NODE_VERSION).tgz: $(NODEFILES)
	npm pack $(DEST)

$(DEST)/geographiclib.js: ../HEADER.js $(SCRIPTS)
	test -d $(DEST) || mkdir $(DEST)
	../js-cat.sh $^ > $@

$(DEST)/geographiclib.min.js: ../HEADER.js $(SCRIPTS)
	test -d $(DEST) || mkdir $(DEST)
	../js-compress.sh $^ > $@

$(DEST)/%.html: ../%.in
	cp $^ $@

$(DEST)/src/%.js: %.js
	test -d $(DEST)/src || mkdir $(DEST)/src
	cp $^ $@

samples: $(DEST)/geod-calc.html $(DEST)/geod-google.html $(DEST)/geod-google-instructions.html

test:
	cd $(DEST) && mocha
