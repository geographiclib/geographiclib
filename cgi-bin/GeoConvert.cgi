#!/bin/sh
# $Id$
. ./utils.sh
OPTION=`lookupkey "$QUERY_STRING" option`
if test "$OPTION" = Reset; then
    INPUT=
else
    INPUT=`lookupkey "$QUERY_STRING" input`
    FORMAT=`lookupkey "$QUERY_STRING" format`
    ZONE=`lookupkey "$QUERY_STRING" zone`
    PREC=`lookupkey "$QUERY_STRING" prec`
fi
test "$FORMAT" || FORMAT=g
test "$ZONE" || ZONE=-2
test "$PREC" || PREC=0
INPUTENC=`encodevalue "$INPUT"`
COMMAND=GeoConvert
EXECDIR=./exec
test $FORMAT = g || COMMAND="$COMMAND -$FORMAT"
case $ZONE in
    -2 ) ;;
    -1 ) COMMAND="$COMMAND -s";;
    * ) COMMAND="$COMMAND -z $ZONE"
esac
test $PREC = 0 || COMMAND="$COMMAND -p $PREC"
if test "$INPUT"; then
    OUTPUT=`echo $INPUT | $EXECDIR/$COMMAND`
    echo `date +"%F %T"` echo $INPUT \| $COMMAND >> ../persistent/utilities.log
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
      Online geographic coordinate conversions
    </title>
  </header>
  <body>
    <h3>
      Online geographic coordinate conversions using the GeoConvert utility
    </h3>
    <form action="/cgi-bin/GeoConvert" method="get">
      <p>
        Location (ex. "<tt>33.3N 44.4E</tt>", "<tt>38SMB4488</tt>", "<tt>38N 444000 3688000</tt>"):<br>
        &nbsp;&nbsp;&nbsp;
        <input type=text name="input" size=30 value="$INPUTENC">
      </p>
      <table>
        <tr>
          <td rowspan="2">
            Output format:<br>
EOF
(
    cat <<EOF
g Decimal degrees
d Degrees minutes seconds
u <a href="http://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system">UTM</a> or <a href="http://en.wikipedia.org/wiki/Universal_Polar_Stereographic">UPS</a>
m <a href="http://en.wikipedia.org/wiki/Military_grid_reference_system">MGRS</a>
EOF
) | while read c desc; do
    CHECKED=
    test "$c" = "$FORMAT" && CHECKED=CHECKED
    echo "&nbsp;&nbsp;&nbsp;"
    echo \<input type=\"radio\" name=\"format\" value=\"$c\" $CHECKED\> $desc\<br\>
done
cat <<EOF
          </td>
          <td>
            &nbsp;&nbsp;&nbsp;
            Output zone:<br>
            &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
            <select name="zone" size=1>
EOF
for ((z=-2; z<=60; ++z)); do
    SELECTED=
    test "$z" = "$ZONE" && SELECTED=SELECTED
    case $z in
        -2 ) name="Match input or standard";;
        -1 ) name="Convert to standard zone";;
        0 ) name=UPS;;
        * ) name="UTM zone $z";;
    esac
    echo \<option $SELECTED value=\"$z\"\>$name
done
cat <<EOF
	    </select>
          </td>
        </tr>
        <tr>
          <td>
            &nbsp;&nbsp;&nbsp;
	    Output precision:<br>
            &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
            <select name="prec" size=1>
EOF
(
    cat <<EOF
-5 100km 1d
-4 10km 0.1d
-3 1km 0.01d 1'
-2 100m 0.001d 0.1'
-1 10m 0.0001d 1"
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
10 0.00000000001"
EOF
) | while read p desc; do
    SELECTED=
    test "$p" = "$PREC" && SELECTED=SELECTED
    echo \<option $SELECTED value=\"$p\"\>$desc
done
cat <<EOF
	    </select>
          </td>
        </tr>
      </table>
      <p>
        Select action:<br>
        &nbsp;&nbsp;&nbsp;
        <input type="submit" name="option" value="Submit">
        <input type="submit" name="option" value="Reset">
      </p>
      <p>
        Results:<br>
        <pre>
    Command = `test "$INPUT" && echo echo $INPUTENC \| $COMMAND`
    Output  = $OUTPUTENC
        </pre>
      </p>
    </form>
    <hr>
    <p>
      <a href="http://geographiclib.sourceforge.net/html/utilities.html#geoconvert">
        GeoConvert</a>
      is one of the utilities provided
      with <a href="http://geographiclib.sourceforge.net/">
        GeographicLib</a>.
      It converts between geographic (latitude and longitude) coordinates,
      <a href="http://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system">
        universal transverse Mercator (UTM)</a> or
      <a href="http://en.wikipedia.org/wiki/Universal_Polar_Stereographic">
        universal polar stereographic (UPS)</a> coordinates, and the
      <a href="http://en.wikipedia.org/wiki/Military_grid_reference_system">
        military grid reference system (MGRS)</a>.
      This web interface illustrates a subset of its capabilities.  If
      you wish to use GeoConvert directly,
      <a href="http://sourceforge.net/projects/geographiclib/files/distrib">
        download</a>
      and compile GeographicLib.  GeoConvert is a simple wrapper in the
      <a href="http://geographiclib.sourceforge.net/html/classGeographicLib_1_1GeoCoords.html">
	GeographicLib::GeoCoords</a> class.
    </p>
    <p>
      Examples of legal geographic locations are (these all refer to the
      same place, Cape Morris Jesup on the northern tip of Greenland):
      <pre>
    Latitude and longitude:       MGRS:
        83.627 -32.664                24XWT783908
        W32d40 N83d37.6               YUB17770380
        83d37'39"N 32d39'52"W     UTM:
    UPS:                              25N 504158 9286521
        N 1617772 1403805             430000 9290000 26N
      </pre>
      The letter in following the zone number in the UTM position is a
      hemisphere designator (N or S) and <em>not</em> the MGRS latitude
      band letter.
    </p>
    <hr>
    <address><a href="http://charles.karney.info/">Charles Karney</a>
      <a href="mailto:charles@karney.com">&lt;charles@karney.com&gt;</a>
      (2009-10-25)</address>
  </body>
</html>
EOF
