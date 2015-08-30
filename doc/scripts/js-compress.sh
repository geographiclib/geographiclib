#! /bin/sh -e
# Concatenate and strip JavaScript files
HEADER=$1
shift
JS_VERSION=`grep -h "version_string = " "$@" | cut -f2 -d'"'`
FILE_INVENTORY=""
for f; do
    FILE_INVENTORY="$FILE_INVENTORY `basename $f`"
done
sed -e "s/@JS_VERSION@/$JS_VERSION/" -e "s/@FILE_INVENTORY@/$FILE_INVENTORY/" \
    $HEADER
for f; do
    echo "// GeographicLib/`basename $f`"
    cat $f | sed -e '1,/\*\//d' -e 's%//.*%%' | tr -s '	 ' ' ' |
    sed -e 's/^ //' -e 's/ $//' | grep -v '^$' |
    sed -e 's/\([^"A-Za-z0-9_]\) /\1/g' -e 's/ \([^\["A-Za-z0-9_]\)/\1/g'
done
