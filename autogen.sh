#!/bin/sh
#
# A simple script to add/reset autotools support
# $Id$

run ()
{
    echo "Running: $*"
    eval $*

    if test $? != 0 ; then
	echo "Error: while running '$*'"
	exit 1
    fi
}

rm -rf autom4ate.cache m4
mkdir m4
run aclocal
run autoheader
run libtoolize --force --copy
run automake --add-missing --copy --force-missing --foreign
run autoconf
run autoreconf

rm -f man/*.usage man/*.1
