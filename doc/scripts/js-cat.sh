#! /bin/sh -e
# Concatenate JavaScript files
HEADER=$1
shift
JS_VERSION=`grep -h "version_string = " "$@" | cut -f2 -d'"'`
FILE_INVENTORY=
for f; do
    FILE_INVENTORY="$FILE_INVENTORY `basename $f`"
done
sed -e "s/@JS_VERSION@/$JS_VERSION/" -e "s/@FILE_INVENTORY@/$FILE_INVENTORY/" \
    $HEADER
for f; do
    echo
    echo "/**************** `basename $f` ****************/"
    cat $f
done
echo
echo "/******** support loading with node's require ********/"
echo "if (typeof module === 'object') module.exports = GeographicLib;"
