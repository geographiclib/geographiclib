VERSION = 1.45
SCRIPTS = Math.js Geodesic.js GeodesicLine.js PolygonArea.js DMS.js
TUTORIALS = tutorials/geodesics.md tutorials/interface.md \
	tutorials/examples.md tutorials/tutorials.json
DEST=node_modules/geographiclib
all: doc script samples
doc: html/index.html
script: $(DEST)/geographiclib.js

# requires that jsdoc be installed
html/index.html: GeographicLib.md $(SCRIPTS) $(TUTORIALS)
	jsdoc -d html -u tutorials -R GeographicLib.md $(SCRIPTS)

publish: doc
	rsync --delete -a -e ssh html/ karney,geographiclib@web.sourceforge.net:htdocs/$(VERSION)/js/

$(DEST)/geographiclib.js: ../HEADER.js $(SCRIPTS)
	test -d $(DEST) || mkdir $(DEST)
	../js-compress.sh $^ > $@

$(DEST)/%.html: ../%.in
	cp $^ $@

samples: $(DEST)/geod-calc.html $(DEST)/geod-google.html $(DEST)/geod-google-instructions.html

test:
	cd $(DEST) && mocha
