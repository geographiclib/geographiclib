#! /bin/sh
# Run Inno Setup Compiler to create installers for geoid datasets

INNO="c:/Program Files/Inno Setup 5/ISCC.exe"
test -f "$INNO" || INNO="c:/Program Files (x86)/Inno Setup 5/ISCC.exe"

GRAVITYDIR=..
test -d "$GRAVITYDIR"/gravity-installers || mkdir -p "$GRAVITYDIR"/gravity-installers
GRAVITYDIR=`cygpath -w $GRAVITYDIR`
(
cat <<EOF
egm84 EGM84
egm96 EGM96
egm2008 EGM2008
wgs84 WGS84
EOF
) | while read prefix name; do
    "$INNO" gravity-installers.iss \
	/dGRAVITYDIR="$GRAVITYDIR" \
	/dPREFIX="$prefix" \
	/dNAME="$name" > $prefix.log 2>&1
done
exit 0
