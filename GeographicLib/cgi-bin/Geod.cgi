#!/bin/sh
# $Id$
. ./utils.sh
OPTION=`lookupkey "$QUERY_STRING" option`
if test "$OPTION" = Reset; then
    INPUT=
else
    INPUT=`lookupcheckkey "$QUERY_STRING" input`
    FORMAT=`lookupkey "$QUERY_STRING" format`
    AZI2=`lookupkey "$QUERY_STRING" azi2`
    PREC=`lookupkey "$QUERY_STRING" prec`
    TYPE=`lookupkey "$QUERY_STRING" type`
fi
test "$FORMAT" || FORMAT=g
test "$AZI2" || AZI2=f
test "$PREC" || PREC=3
test "$TYPE" || TYPE=d
INPUTENC=`encodevalue "$INPUT"`
COMMAND=Geod
EXECDIR=./exec
test $TYPE = d || COMMAND="$COMMAND -$TYPE"
test $FORMAT = g || COMMAND="$COMMAND -$FORMAT"
test $AZI2 = f || COMMAND="$COMMAND -$AZI2"
test $PREC = 3 || COMMAND="$COMMAND -p $PREC"
if test "$INPUT"; then
    OUTPUT=`echo $INPUT | $EXECDIR/$COMMAND`
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
      Online geodesic calculations
    </title>
  </header>
  <body>
    <h3>
      Online geodesic calculations using the Geod utility
    </h3>
    <form action="/cgi-bin/Geod" method="get">
      <p>
        Input (ex. "<tt>40N 74W 53.2 5900e3</tt>" [direct], "<tt>40N 74W 49N 2E</tt>" [inverse]):<br>
        &nbsp;&nbsp;&nbsp;
        <input type=text name="input" size=50 value="$INPUTENC">
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
    test "$c" = "$AZI2" && CHECKED=CHECKED
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
10 100pm 0.00000000001"
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
        Geodesic calculation:
	<table>
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
	      &rarr; <em>lat2 lon2 azi2 m12</em>
	    </td>
	  <tr>
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
	      &rarr; <em>azi1 azi2 s12 m12</em>
	    </td>
	  <tr>
	</table>
      </p>
      <p>
        Select action:<br>
        &nbsp;&nbsp;&nbsp;
        <input type="submit" name="option" value="Submit">
        <input type="submit" name="option" value="Reset">
      </p>
      <p>
        Results:<br>
        <pre>
    Command = `test "$INPUT" && echo "echo $INPUTENC | $COMMAND"`
    Output  = $OUTPUTENC
        </pre>
      </p>
    </form>
    <hr>
    <p>
      <a href="http://geographiclib.sourceforge.net/html/utilities.html#geod">
        Geod</a>
      is one of the utilities provided
      with <a href="http://geographiclib.sourceforge.net/">
        GeographicLib</a>.
      It performs direct and inverse geodesic calculations for the WGS84
      ellipsoid.  The shortest path between two points on the ellipsoid
      at (<em>lat1</em>, <em>lon1</em>) and (<em>lat2</em>,
      <em>lon2</em>) is called the geodesic.  Its length is <em>s12</em>
      and the geodesic from point 1 to point 2 has azimuths
      <em>azi1</em> and <em>azi2</em> at the two end points.  The
      reduced length of the geodesic, <em>m12</em>, is defined such that
      if the initial azimuth is perturbed by <em>dazi1</em> (radians)
      then the second point is displaced by <em>m12</em> *
      <em>dazi1</em> in the direction perpendicular to the geodesic.
      There are two standard geodesic problems:
      <ul>
	<li> Direct: the input is <em>lat1 lon1 azi1 s12</em> and the
	  output is <em>lat2 lon2 azi2 m12</em>;
	<li> Inverse: the input is <em>lat1 lon1 lat2 lon2</em> and the
	  output is <em>azi1 azi2 s12 m12</em>.
      </ul>
      This web interface illustrates a subset of the capabilities of
      Geod.  If you wish to use Geod directly,
      <a href="http://sourceforge.net/projects/geographiclib/files/distrib">
        download</a>
      and compile GeographicLib.  Geod is a simple wrapper in the
      <a href="http://geographiclib.sourceforge.net/html/classGeographicLib_1_1Geodesic.html">
	GeographicLib::Geodesic</a> class.
      A description of the algorithms is given
      <a href="http://geographiclib.sourceforge.net/html/geodesic.html">
	here</a>.
    </p>
    <p>
      Latitudes and longitudes can be given in various formats, for
      example (these all refer to the position of Timbuktu):
      <pre>
        16.776 -3.009
        16d47' -3d1'
        W3d0'34" N16d46'33"
      </pre>
      Azimuths are given in degress clockwise from north.  The
      distances <em>s12</em> and <em>m12</em> are in meters.
    </p>
    <p>
      The "standard" method of calculating geodesics uses an algorithm
      given by
      <a href="http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf">
	Vincenty (1975)</a>.
      However, this has limited accuracy and completely fails for the
      inverse problem when the points are nearly antipodal.  The method
      used by Geod is accurate to about 15 nm and gives solutions for
      the inverse problem for any pair of points.  For example, compare
      the inverse result given by Geod for the antipodal points (N30,
      E0) and (S30, E180) where the geodesic follows a meridian with the
      bogus result returned by the
      <a href="http://www.ngs.noaa.gov/">
	NGS</a> online 
      <a href="http://www.ngs.noaa.gov/cgi-bin/Inv_Fwd/inverse2.prl">
	inverse geodesic calculator</a>.
    </p>
    <hr>
    <address><a href="http://charles.karney.info/">Charles Karney</a>
      <a href="mailto:charles@karney.com">&lt;charles@karney.com&gt;</a>
      (2009-10-26)</address>
  </body>
</html>
EOF
