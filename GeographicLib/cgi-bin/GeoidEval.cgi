#!/bin/sh
# $Id$
. ./utils.sh
OPTION=`lookupkey "$QUERY_STRING" option`
if test "$OPTION" = Reset; then
    INPUT=
else
    INPUT=`lookupcheckkey "$QUERY_STRING" input`
    GEOID=`lookupkey "$QUERY_STRING" geoid`
fi
test "$GEOID" || GEOID=egm96-5
INPUTENC=`encodevalue "$INPUT"`
COMMAND=GeoidEval
GEOID_PATH=geoids
EXECDIR=./exec
test $GEOID = egm96-5 || COMMAND="$COMMAND -n $GEOID"
if test "$INPUT"; then
    OUTPUT=`echo $INPUT | GEOID_PATH=$GEOID_PATH $EXECDIR/$COMMAND`
    test $? -eq 0 && OUTPUT="`echo $OUTPUT | cut -f1 -d' '` m"
    echo `date +"%F %T"` "echo $INPUT | $COMMAND" >> ../persistent/utilities.log
else
    OUTPUT=
    echo `date +"%F %T"` $COMMAND >> ../persistent/utilities.log
fi
OUTPUTENC=`encodevalue "$OUTPUT"`

echo Content-type: text/html
echo
cat <<EOF
<html>
  <header>
    <title>
      Online geoid calculator
    </title>
  </header>
  <body>
    <h3>
      Online geoid calculator using the GeoidEval utility
    </h3>
    <form action="/cgi-bin/GeoidEval" method="get">
      <p>
	Position (ex. "<tt>16.8N 3.0W</tt>"):<br>
	&nbsp;&nbsp;&nbsp;
	<input type=text name="input" size=30 value="$INPUTENC">
      </p>
      <p>
	Earth gravity model:<br>
EOF
(
    cat <<EOF
egm84-15 <a href="http://earth-info.nga.mil/GandG/wgs84/gravitymod/wgs84_180/wgs84_180.html">EGM84</a>
egm96-5 <a href="http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html">EGM96</a>
egm2008-1 <a href="http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008">EGM2008</a>
EOF
) | while read c desc; do
    CHECKED=
    test "$c" = "$GEOID" && CHECKED=CHECKED
    echo "&nbsp;&nbsp;&nbsp;"
    echo "<input type=\"radio\" name=\"geoid\" value=\"$c\" $CHECKED> $desc"
done
cat <<EOF
      </p>
      <p>
	Select action:<br>
	&nbsp;&nbsp;&nbsp;
	<input type="submit" name="option" value="Submit">
	<input type="submit" name="option" value="Reset">
      </p>
      <p>
	Geoid height:<br>
	<pre>
    Command = `test "$INPUT" && echo "echo $INPUTENC | $COMMAND"`
    Height  = $OUTPUTENC
	</pre>
      </p>
    </form>
    <hr>
    <p>
      <a href="http://geographiclib.sourceforge.net/html/utilities.html#geoideval">
	GeoidEval</a>
      is one of the utilities provided
      with <a href="http://geographiclib.sourceforge.net/">
	GeographicLib</a>.
      It computes the height of the geoid above the WGS84 ellipsoid
      using interpolation in a grid of values for one of the earth
      gravity models:
      <ul>
	<li>
	  <a href="http://earth-info.nga.mil/GandG/wgs84/gravitymod/wgs84_180/wgs84_180.html">
	    EGM84</a>
	<li>
	  <a href="http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html">
	    EGM96</a>
	<li>
	  <a href="http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008">
	    EGM2008</a>
      </ul>
      The RMS error in the interpolated height is about 1 mm.  This web
      interface illustrates a subset of its capabilities.  If you wish
      to use GeoidEval directly,
      <a href="http://sourceforge.net/projects/geographiclib/files/distrib">
	download</a>
      and compile GeographicLib.  GeoidEval is a simple wrapper in the
      <a href="http://geographiclib.sourceforge.net/html/classGeographicLib_1_1Geoid.html">
	GeographicLib::Geoid</a> class.
      A description of the methods is given
      <a href="http://geographiclib.sourceforge.net/html/geoid.html">
	here</a>.
    </p>
    <p>
      Give the position in terms of latitude and longitude, for example
      (these all refer to the position of Timbuktu):
      <pre>
        16.776 -3.009
        16d47' -3d1'
        W3d0'34" N16d46'33"
      </pre>
    </p>
    <hr>
    <address><a href="http://charles.karney.info/">Charles Karney</a>
      <a href="mailto:charles@karney.com">&lt;charles@karney.com&gt;</a>
      (2009-10-25)</address>
  </body>
</html>
EOF
