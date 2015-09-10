SCRIPTS = Math.js Geodesic.js GeodesicLine.js PolygonArea.js DMS.js
TUTORIALS = tutorials/geodesics.md tutorials/interface.md tutorials/examples.md
DEST=node_modules/geographiclib
all: doc script samples
doc: html/index.html
script: $(DEST)/geographiclib.js

# requires that jsdoc and ink-docstrap are installed
html/index.html: README.md $(SCRIPTS) $(TUTORIALS)
	jsdoc -d html -t /usr/lib/node_modules/ink-docstrap/template -c conf.json -u tutorials README.md Math.js Geodesic.js GeodesicLine.js PolygonArea.js DMS.js

$(DEST)/geographiclib.js: ../HEADER.js $(SCRIPTS)
	test -d $(DEST) || mkdir $(DEST)
	../js-compress.sh $^ > $@

$(DEST)/%.html: ../%.in
	cp $^ $@

samples: $(DEST)/geod-calc.html $(DEST)/geod-google.html $(DEST)/geod-google-instructions.html

test:
	cd $(DEST) && mocha
