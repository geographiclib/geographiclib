SCRIPTS = Math.js Geodesic.js GeodesicLine.js PolygonArea.js DMS.js
TUTORIALS = tutorials/geodesics.md tutorials/interface.md tutorials/examples.md tutorials/tutorials.json
DEST=node_modules/geographiclib
all: doc script samples
doc: html/index.html
script: $(DEST)/geographiclib.js

# requires that jsdoc and ink-docstrap are installed
html/index.html: GeographicLib.md $(SCRIPTS) $(TUTORIALS)
	jsdoc -d html -t /usr/lib/node_modules/ink-docstrap/template -c conf.json -u tutorials -R GeographicLib.md $(SCRIPTS)

publish: doc
#	rsync --delete -av -e ssh html/ vorpal2:/var/tmp/html/
	rsync --delete -av -e ssh html/ karney,geographiclib@web.sourceforge.net:htdocs/1.45/js/

$(DEST)/geographiclib.js: ../HEADER.js $(SCRIPTS)
	test -d $(DEST) || mkdir $(DEST)
	../js-compress.sh $^ > $@

$(DEST)/%.html: ../%.in
	cp $^ $@

samples: $(DEST)/geod-calc.html $(DEST)/geod-google.html $(DEST)/geod-google-instructions.html

test:
	cd $(DEST) && mocha
