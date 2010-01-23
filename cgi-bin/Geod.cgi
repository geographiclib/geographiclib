#!/bin/sh
# $Id$
. ./utils.sh
OPTION=`lookupkey "$QUERY_STRING" option`
if test "$OPTION" = Reset; then
    INPUT=
else
    INPUT=`lookupcheckkey "$QUERY_STRING" input`
    FORMAT=`lookupkey "$QUERY_STRING" format`
    AZF2=`lookupkey "$QUERY_STRING" azi2`
    PREC=`lookupkey "$QUERY_STRING" prec`
    TYPE=`lookupkey "$QUERY_STRING" type`
fi
test "$FORMAT" || FORMAT=g
test "$AZF2" || AZF2=f
test "$PREC" || PREC=3
test "$TYPE" || TYPE=i
AZX="faz2"
test "$AZF2" = b && AZX="baz2"

INPUTENC=`encodevalue "$INPUT"`
COMMAND=Geod
EXECDIR=../bin
F='<font color="blue">'
G='</font>'
test $TYPE = d || COMMAND="$COMMAND -$TYPE"
COMMANDX="$COMMAND -f -p 1"
test $FORMAT = g || COMMAND="$COMMAND -$FORMAT"
test $AZF2 = f || COMMAND="$COMMAND -$AZF2"
test $PREC = 3 || COMMAND="$COMMAND -p $PREC"
STATUS=
POSITION1=
POSITION2=
DIST12=
if test "$INPUT"; then
    COMMANDLINE="echo $INPUT | $COMMAND"
    OUTPUT=`echo $INPUT | $EXECDIR/$COMMAND -f`
    if test $? -eq 0; then
	STATUS=OK
	OUTPUTG=`echo $INPUT | $EXECDIR/$COMMANDX`
	POS1="`echo $OUTPUT | cut -f1-2 -d' '`"
	POS2="`echo $OUTPUT | cut -f4-5 -d' '`"
	POSG1="`echo $OUTPUTG | cut -f1-2 -d' '`"
	POSG2="`echo $OUTPUTG | cut -f4-5 -d' '`"
	AZI1="`echo $OUTPUT | cut -f3 -d' '`"
	AZI2="`echo $OUTPUT | cut -f6 -d' '`"
	DIST12="`echo $OUTPUT | cut -f7 -d' '`"
	if test "$TYPE" = d; then
	    POSITION1=$(geohack $POSG1 $POS1 Black)\ $(encodevalue "$AZI1")
	    POSITION2=$F$(geohack $POSG2 $POS2 Blue)\ $(encodevalue "$AZI2")$G
	    DIST12=$(encodevalue "$DIST12")
	else
	    POSITION1=$(geohack $POSG1 $POS1 Black)\ $F$(encodevalue "$AZI1")$G
	    POSITION2=$(geohack $POSG2 $POS2 Black)\ $F$(encodevalue "$AZI2")$G
	    DIST12=$F$(encodevalue "$DIST12")$G
	fi
    else
	STATUS="$OUTPUT"
    fi
    # echo `date +"%F %T"` "$COMMAND: $INPUT" >> ../persistent/utilities.log
fi

echo Content-type: text/html
echo
cat <<EOF
<html>
  <head>
    <title>
      Online geodesic calculator
    </title>
    <meta name="description" content="Online geodesic calculator" />
    <meta name="author" content="Charles F. F. Karney" />
    <meta name="keywords"
	  content="geodesics,
		   geodesic distance,
		   geographic distance,
		   shortest path,
		   direct geodesic problem,
		   inverse geodesic problem,
		   distance and azimuth,
		   distance and heading,
		   range and bearing,
		   spheroidal triangle,
		   latitude and longitude,
		   online calculator,
		   WGS84 ellipsoid,
		   GeographicLib" /> 
  </head>
  <body>
    <h3>
      Online geodesic calculations using the
      <a href="http://geographiclib.sourceforge.net/html/utilities.html#geod">
	 Geod</a> utility
    </h3>
    <form action="/cgi-bin/Geod" method="get">
      <p>
        Geodesic calculation:
        <table>
          <tr>
            <td>
              &nbsp;&nbsp;&nbsp;
              <input type="radio" name="type" value="i"
                     `test "$TYPE" = i && echo CHECKED`>
            </td>
            <td>
              Inverse:
            </td>
            <td>
              <em>lat1 lon1 lat2 lon2</em>
            </td>
            <td>
              &rarr; <em>azi1 azi2 s12</em>
            </td>
          <tr>
          <tr>
            <td>
              &nbsp;&nbsp;&nbsp;
              <input type="radio" name="type" value="d"
                     `test "$TYPE" = d && echo CHECKED`>
            </td>
            <td>
              Direct:
            </td>
            <td>
              <em>lat1 lon1 azi1 s12</em>
            </td>
            <td>
              &rarr; <em>lat2 lon2 azi2</em>
            </td>
          <tr>
        </table>
      </p>
      <p>
        Input (ex. &laquo;<tt>40.6 -73.8 49d01'N 2d33'E</tt>&raquo; [inverse],
	&laquo;<tt>40d38'23"N 073d46'44"W 53d30' 5850e3</tt>&raquo; [direct]):
	<br>
        &nbsp;&nbsp;&nbsp;
        <input type=text name="input" size=72 value="$INPUTENC">
      </p>
      <p>
        <table>
          <tr>
            <td>
              Output format:
            </td>
EOF
(
    cat <<EOF
g Decimal degrees
d Degrees minutes seconds
EOF
) | while read c desc; do
    CHECKED=
    test "$c" = "$FORMAT" && CHECKED=CHECKED
    echo "<td>&nbsp;"
    echo "<input type=\"radio\" name=\"format\" value=\"$c\" $CHECKED> $desc"
    echo "</td>"
done
cat <<EOF
          </tr>
          <tr>
            <td>
              Heading at point 2:
            </td>
EOF
(
    cat <<EOF
f Forward azimuth
b Back azimuth
EOF
) | while read c desc; do
    CHECKED=
    test "$c" = "$AZF2" && CHECKED=CHECKED
    echo "<td>&nbsp;"
    echo "<input type=\"radio\" name=\"azi2\" value=\"$c\" $CHECKED> $desc"
    echo "</td>"
done
cat <<EOF
          </tr>
          <tr>
            <td>
              Output precision:
            </td>
            <td colspan="2">&nbsp;
              <select name="prec" size=1>
EOF
(
    cat <<EOF
0 1m 0.00001d 0.1"
1 100mm 0.01"
2 10mm 0.001"
3 1mm 0.0001"
4 100um 0.00001"
5 10um 0.000001"
6 1um 0.0000001"
7 100nm 0.00000001"
8 10nm 0.000000001"
9 1nm 0.0000000001"
EOF
) | while read p desc; do
    SELECTED=
    test "$p" = "$PREC" && SELECTED=SELECTED
    echo "<option $SELECTED value=\"$p\"> $desc<br>"
done
cat <<EOF
              </select>
            </td>
          </tr>
        </table>
      </p>
      <p>
        Select action:<br>
        &nbsp;&nbsp;&nbsp;
        <input type="submit" name="option" value="Submit">
        <input type="submit" name="option" value="Reset">
      </p>
      <p>
        Geodesic (input in black, output in ${F}blue${G}):<br>
        <font size="4"><pre>
    status         = `encodevalue "$STATUS"`
    lat1 lon1 faz1 = $POSITION1
    lat2 lon2 $AZX = $POSITION2
    s12 (m)        = $DIST12</pre></font>
      </p>
    </form>
    <hr>
    <p>
      <a href="http://geographiclib.sourceforge.net/html/utilities.html#geod">
        Geod</a>
      performs geodesic calculations for the WGS84 ellipsoid.  The
      shortest path between two points on the ellipsoid at
      (<em>lat1</em>, <em>lon1</em>) and (<em>lat2</em>,
      <em>lon2</em>) is called the geodesic; its length is <em>s12</em>
      and the geodesic from point 1 to point 2 has azimuths
      <em>azi1</em> and <em>azi2</em> at the two end points.
      There are two standard geodesic problems:
      <ul>
        <li> Direct: &nbsp; given [<em>lat1 lon1 azi1 s12</em>] find
          [<em>lat2 lon2 azi2</em>];
        <li> Inverse: given [<em>lat1 lon1 lat2 lon2</em>] find
          [<em>azi1 azi2 s12</em>].
      </ul>
      Latitudes and longitudes can be given in various formats, for
      example (these all refer to the position of Timbuktu):
      <pre>
        16.776 -3.009
        16d47' -3d1'
        W3d0'34" N16d46'33"</pre>
      Azimuths are given in degress clockwise from north.  The
      distance <em>s12</em> is in meters.
    </p>
    <p>
      The "standard" method of calculating geodesics uses an algorithm
      given by
      <a href="http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf">
        Vincenty (1975)</a>.
      However, this has limited accuracy and fails for the inverse
      problem when the points are nearly antipodal.  The method used by
      Geod is accurate to about 15 nm and gives solutions for the
      inverse problem for any pair of points.  For example, compare the
      inverse result given by Geod for the antipodal points (N30, E0)
      and (S30, E180) where the geodesic follows a meridian with the
      bogus result returned by the
      <a href="http://www.ngs.noaa.gov/">
        NGS</a> online 
      <a href="http://www.ngs.noaa.gov/cgi-bin/Inv_Fwd/inverse2.prl">
        inverse geodesic calculator</a>.
    </p>
    <p>
      <a href="http://geographiclib.sourceforge.net/html/utilities.html#geod">
        Geod</a>,
      which is a simple wrapper of the
      <a href="http://geographiclib.sourceforge.net/html/classGeographicLib_1_1Geodesic.html">
        GeographicLib::Geodesic</a> class,
      is one of the utilities provided
      with <a href="http://geographiclib.sourceforge.net/">
        GeographicLib</a>.
      This web interface illustrates a subset of its capabilities.  If
      you wish to use Geod directly,
      <a href="http://sourceforge.net/projects/geographiclib/files/distrib">
        download</a>
      and compile GeographicLib.  A description of the algorithms is
      given
      <a href="http://geographiclib.sourceforge.net/html/geodesic.html">
        here</a>.
    </p>
    <hr>
    <address><a href="http://charles.karney.info/">Charles Karney</a>
      <a href="mailto:charles@karney.com">&lt;charles@karney.com&gt;</a>
      (2010-01-08)</address>
    <a href="http://sourceforge.net">
      <img src="http://sourceforge.net/sflogo.php?group_id=283628&amp;type=2" border="0" alt="SourceForge.net" />
    </a>
  </body>
</html>
EOF
