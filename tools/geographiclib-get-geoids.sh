#! /bin/sh
#
# Downloads geoid datasets for use by GeographicLib::Geoid.  This is
# modeled on a similar script geographiclib-datasets-download by
# Francesco P. Lovergine <frankie@debian.org>
#
# Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed
# under the MIT/X11 License.  For more information, see
# http://geographiclib.sourceforge.net/
#
# $Id$

GEOID_DEFAULT_PATH="@GEOID_DEFAULT_PATH@"
DEFAULTDIR=`dirname "$GEOID_DEFAULT_PATH"`
usage() {
    cat <<EOF
usage: $0 [-p parentdir] [-d] [-h] geoid...

This program downloads and installs the datasets used by the
GeographicLib::Geoid class and the GeoidEval tool to compute geoid
heights.  These datasets are NGA earth gravity models evaluated on a
rectangular grid in latitude and longitude.  geoid is one of more of the
names from this table:

                                  size (MB)  
  name         geoid    grid    tar.bz2  disk
  egm84-30     EGM84    30'      0.5      0.6
  egm84-15     EGM84    15'      1.5      2.1
  egm96-15     EGM96    15'      1.5      2.1
  egm96-5      EGM96     5'       11       19
  egm2008-5    EGM2008   5'       11       19
  egm2008-2_5  EGM2008   2.5'     35       75
  egm2008-1    EGM2008   1'      170      470

The size columns give the download and installed sizes of the datasets.
In addition you can specify

  all = all of the datasets
  minimal = emg96-5
  best = egm84-15 egm96-5 egm2008-1 (the highest resolution for each
         earth gravity model)
  good = same as best but substitute egm2008-2_5 for egm2008-1 to save
         on disk space

If no name is specified then minimal is assumed.

-p parentdir (default $DEFAULTDIR) specifies where the
datasets should be stored.  The "Default geoid path" listed when running

  GeoidEval -h

should be parentdir/geoids.  This script must
be run by a user with write access to this directory.

If -d is provided, the temporary directory which holds the downloads,
${TMPDIR:-/tmp}/geoid-XXXXXXXX, will be saved.  -h prints this help.

For more information on the geoid datasets, visit

  http://geographiclib.sourceforge.net/html/geoid.html

EOF
}

PARENTDIR="$DEFAULTDIR"
DEBUG=
while getopts hp:d c; do
    case $c in
        h )
            usage;
            exit 0
            ;;
        p ) PARENTDIR="$OPTARG"
            ;;
        d ) DEBUG=y
            ;;
        * )
            usage 1>&2;
            exit 1
            ;;
    esac
done
shift `expr $OPTIND - 1`

test -d "$PARENTDIR"/geoids || mkdir -p "$PARENTDIR"/geoids 2> /dev/null
if test ! -d "$PARENTDIR"/geoids; then
    echo Cannot create directory $PARENTDIR/geoids 1>&2
    exit 1
fi

TEMP=
if test -z "$DEBUG"; then
trap 'trap "" 0; test "$TEMP" && rm -rf "$TEMP"; exit 1' 1 2 3 9 15
trap            'test "$TEMP" && rm -rf "$TEMP"'            0
fi
TEMP=`mktemp --tmpdir --quiet --directory geoid-XXXXXXXX`

if test -z "$TEMP" -o ! -d "$TEMP"; then
    echo Cannot create temporary directory 1>&2
    exit 1
fi

WRITETEST="$PARENTDIR"/geoids/write-test-`basename $TEMP`
if touch "$WRITETEST" 2> /dev/null; then
    rm -f "$WRITETEST"
else
    echo Cannot write in directory $PARENTDIR/geoids 1>&2
    exit 1
fi
    
set -e

cat > $TEMP/all <<EOF
egm84-30
egm84-15
egm96-15
egm96-5
egm2008-5
egm2008-2_5
egm2008-1
EOF

test $# -eq 0 && set -- minimal

while test $# -gt 0; do
    if grep "^$1\$" $TEMP/all > /dev/null; then
	echo $1
    else
	case "$1" in
	    all )
		cat $TEMP/all
		;;
	    minimal )		# same as no argument
		echo egm96-5
		;;
	    best )		# highest resolution models
		cat <<EOF
egm2008-1
egm96-5
egm84-15
EOF
		;;
	    good )	   # like best but with egm2008-1 -> egm2008-2_5
		cat <<EOF
egm2008-2_5
egm96-5
egm84-15
EOF
		;;
	    * )
		echo Unknown geoid $1 1>&2
		exit 1
		;;
	esac
    fi
    shift
done > $TEMP/list

sort -u $TEMP/list > $TEMP/todo

while read file; do
    echo download $file.tar.bz2 ...
    URL="http://downloads.sourceforge.net/project/geographiclib/geoids-distrib/$file.tar.bz2?use_mirror=autoselect"
    ARCHIVE=$TEMP/$file.tar.bz2
    wget -O$ARCHIVE $URL
    echo unpack $file.tar.bz2 ...
    tar vxojf $ARCHIVE -C $PARENTDIR
    echo geoid $file installed.
done < $TEMP/todo

if test "$DEBUG"; then
    echo Saving temporary directory $TEMP
fi
cat <<EOF

Geoid datasets `tr '\n' ' ' < $TEMP/todo`
downloaded and installed in $PARENTDIR/geoids.

EOF
