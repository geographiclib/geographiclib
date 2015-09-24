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
cat <<EOF

(function() {
  'use strict';
EOF
for f; do
    echo
    echo "/**************** `basename $f` ****************/"
    cat $f
done
cat <<EOF

  /******** export GeographicLib ********/
  if (typeof module === 'object' && module.exports)
    module.exports = GeographicLib;
  else if (typeof define === 'function' && define.amd)
    define('geographiclib', [], function() { return GeographicLib; });
  else if (typeof window === 'object')
    window.GeographicLib = GeographicLib;
  else
    return GeographicLib;
})();
EOF
