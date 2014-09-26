#! /bin/sh
# Concatenate and strip JavaScript files
cat <<EOF
/*
 * Geodesic routines from GeographicLib translated to JavaScript.  For
 * more information, see
 * http://geographiclib.sf.net/html/other.html#javascript
 *
 * The algorithms are derived in
 *
 *    Charles F. F. Karney,
 *    Algorithms for geodesics, J. Geodesy 87, 43-55 (2013),
 *    http://dx.doi.org/10.1007/s00190-012-0578-z
 *    Addenda: http://geographiclib.sf.net/geod-addenda.html
 *
 * Copyright (c) Charles Karney (2011-2014) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sf.net/
 *
 * Inventory of files;
EOF
for f; do
    echo " *  " `basename $f`
done
echo " */"
for f; do
    echo "// `basename $f`"
    cat $f | sed -e '1,/\*\//d' -e 's%//.*%%' | tr -s '	 ' ' ' |
    sed -e 's/^ //' -e 's/ $//' | grep -v '^$' |
    sed -e 's/\([^",:A-Za-z0-9_]\) /\1/g' |
    sed -e 's/ \([^\[":A-Za-z0-9_]\)/\1/g' |
    sed -e 's/^\([^"]*\) :/\1:/' -e 's/^\([^"]*\) :/\1:/' |
    sed -e 's/^\([^"]*\) :/\1:/' -e 's/^\([^"]*\) :/\1:/' |
    sed -e 's/^\([^"]*\) :/\1:/' -e 's/^\([^"]*\) :/\1:/' |
    sed -e 's/\([:,]\) \([^"]*\)$/\1\2/' -e 's/\([:,]\) \([^"]*\)$/\1\2/' |
    sed -e 's/\([:,]\) \([^"]*\)$/\1\2/' -e 's/\([:,]\) \([^"]*\)$/\1\2/' |
    sed -e 's/\([:,]\) \([^"]*\)$/\1\2/' -e 's/\([:,]\) \([^"]*\)$/\1\2/' |
    sed -e 's/, "/,"/g' -e 's/" : "/":"/g' -e 's/\([^ ]\): /\1:/g'
done
