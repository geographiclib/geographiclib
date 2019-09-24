#! /bin/sh
# Convert geodesic-biblio.txt to an html page
cat <<'EOF'
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN"
          "http://www.w3.org/TR/html4/strict.dtd">
<html>
  <head>
    <title>
      Geodesic bibliography
    </title>
    <link rel="stylesheet" type="text/css" href="../default.css">
    <meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1">
  </head>
  <body topmargin=10 leftmargin=10>
    <h3>A geodesic bibliography</h3>
    <p>
      Here is a list of the older mathematical treatments of the
      geodesic problem for an ellipsoid, together with links to online
      copies.  This includes some more recent works which are available
      online (chiefly works funded by the US Government).  Where
      available, links to translations, especially into English, have
      been added.  Unfortunately, the fold-out pages of figures in some
      books are usually not scanned properly by Google; in some cases I
      have been able to scan the missing pages.  In addition, readers
      may not have access to the full text of some Google Books; in
      those cases, I have provided a "pdf" link.  See the files in <a
      href=".">geodesic-papers</a>. Please let me, Charles Karney
      <a href="mailto:charles@karney.com">&lt;charles@karney.com&gt;</a>,
      know of errors, omissions, etc.  In particular, I'm interested to
      learn of any cases where I have mis-translated the title of a
      paper.  In addition, the links to Google Books occasionally get
      out of date; let me know if this happens.
    </p>
    <p>This bibliography was started on 2009-06-06 (at
      <a href="http://trac.osgeo.org/proj/wiki/GeodesicCalculations">
        http://trac.osgeo.org/proj/wiki/GeodesicCalculations</a>, now
      defunct) and
      moved to this site on 2011-02-01.  The last update was on
EOF
git log --date=short $1 | head -3 | tail -1 | tr -s ' ' '	' |
cut -f2 | sed 's/$/./'
cat <<EOF
    </p>
    <ul>
EOF
cat $1 |
sed -e 's/\*/<li>/' -e 's/  *\[\[BR\]\]/ <br>/' -e 's/\[\[BR\]\]/<br>/' \
    -e "s%'''\([0-9][0-9]*\)'''%<b>\1</b>%g" \
    -e "s% ''% <i>%g" -e "s%\([^ ]\)''%\1</i>%g" \
    -e 's%\(https\?\)://\([a-zA-Z][^ ]*\)%<a href="\1://FIX1\2FIX2">\1://\2</a>%' \
    -e 's/FIX1\([a-zA-Z][^ ]*\)"\([^ ]*\)FIX2/FIX1\1%22\2FIX2/g' \
    -e 's/FIX1\([a-zA-Z][^ ]*\)"\([^ ]*\)FIX2/FIX1\1%22\2FIX2/g' \
    -e 's/FIX[12]//g' \
    -e 's%(PDF \([^)]*\))%(<a href="https://geographiclib.sourceforge.io/geodesic-papers/\1">pdf</a>)%' \
    -e 's/&/\&amp;/g' \
    -e 's/\([0-9]\)--\([0-9]\)/\1\&ndash;\2/g' \
    -e 's/É/\&Eacute;/g' \
    -e 's/é/\&eacute;/g' \
    -e 's/á/\&aacute;/g' \
    -e 's/à/\&agrave;/g' \
    -e 's/è/\&egrave;/g' \
    -e 's/ê/\&ecirc;/g' \
    -e 's/ù/\&ugrave;/g' \
    -e 's/ç/\&ccedil;/g' \
    -e 's/ä/\&auml;/g' \
    -e 's/Ü/\&Uuml;/g' \
    -e 's/ï/\&iuml;/g' \
    -e 's/í/\&iacute;/g' \
    -e 's/ö/\&ouml;/g' \
    -e 's/ü/\&uuml;/g' \
    -e 's/ß/\&szlig;/g' | awk '
BEGIN {
    quote=0;
}
{
    if ($0 ~ /^    /) {
        if (!quote) {
            printf "<blockquote>";
            quote = 1;
        }
    } else {
        if (quote) {
            printf "</blockquote>";
            quote = 0;
        }
    }
    print $0;
}
END {
    if (quote)
        printf "</p></blockquote>/n";
}'
cat <<'EOF'
    </ul>
    <hr>
    <a href="..">GeographicLib home</a>
  </body>
</html>
EOF
