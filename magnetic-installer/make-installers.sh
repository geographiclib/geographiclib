#! /bin/sh
# Run Inno Setup Compiler to create installers for geoid datasets
# $Id$
INNO="c:/Program Files/Inno Setup 5/ISCC.exe"
test -f "$INNO" || INNO="c:/Program Files (x86)/Inno Setup 5/ISCC.exe"

MAGNETICDIR=..
test -d "$MAGNETICDIR"/magnetic-installers || mkdir -p "$MAGNETICDIR"/magnetic-installers
MAGNETICDIR=`cygpath -w $MAGNETICDIR`
(
cat <<EOF
wmm2010 WMM2010
emm2010 EMM2010
igrf11  IGRF11
EOF
) | while read prefix name; do
    "$INNO" magnetic-installers.iss \
	/dMAGNETICDIR="$MAGNETICDIR" \
	/dPREFIX="$prefix" \
	/dNAME="$name" > $prefix.log 2>&1
done
exit 0
