#! /bin/sh
# Concatenate and strip Javascript files
cat <<EOF
/*
 * Geodesic routines from GeographicLib translated to Javascript.  For
 * more information, see
 * http://geographiclib.sf.net/html/other.html#javascript
 *
 * Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sf.net/
 *
 * Inventory of files;
EOF
for f; do
    id=`ident < $f | sed -e 's/^ *//'`
    echo " *  " `basename $f` $id
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

