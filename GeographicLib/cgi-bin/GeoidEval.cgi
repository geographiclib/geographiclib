#!/bin/sh
# $Id$
. ./utils.sh
OPTION=`lookupkey "$QUERY_STRING" option`
if test "$OPTION" = Reset; then
    INPUT=
else
    INPUT=`lookupcheckkey "$QUERY_STRING" input`
fi
INPUTENC=`encodevalue "$INPUT"`
COMMAND=GeoidEval
GEOID_PATH=../geoids
EXECDIR=../bin
F='<font color="blue">'
G='</font>'
POSITION1=
POSITION2=
HEIGHT96=
HEIGHT84=
HEIGHT2008=
if test "$INPUT"; then
    HEIGHT96=`echo $INPUT |
    GEOID_PATH=$GEOID_PATH $EXECDIR/$COMMAND -n egm96-5`
    if test $? -eq 0; then
	POSITION1=`echo $INPUT | $EXECDIR/GeoConvert`
	POSITION1=`geohack $POSITION1 $POSITION1 Black`
	POSITION2=\(`echo $INPUT | $EXECDIR/GeoConvert -d -p -1`\)
	HEIGHT2008=`echo $INPUT |
	GEOID_PATH=$GEOID_PATH $EXECDIR/$COMMAND -n egm2008-1`
	HEIGHT84=`echo $INPUT |
	GEOID_PATH=$GEOID_PATH $EXECDIR/$COMMAND -n egm84-15`
	HEIGHT2008=`echo $HEIGHT2008 | cut -f1 -d' '`
	HEIGHT96=`echo $HEIGHT96 | cut -f1 -d' '`
	HEIGHT84=`echo $HEIGHT84 | cut -f1 -d' '`
    else
	POSITION1=`encodevalue "$HEIGHT96"`
	HEIGHT96=
    fi
    # echo `date +"%F %T"` "$COMMAND: $INPUT" >> ../persistent/utilities.log
fi

echo Content-type: text/html
echo
cat <<EOF
<html>
  <head>
    <title>
      Online geoid calculator
    </title>
    <meta name="description" content="Online geoid calculator" />
    <meta name="author" content="Charles F. F. Karney" />
    <meta name="keywords"
	  content="geoid height,
		   orthometric height,
		   earth gravity model,
		   EGM84, EGM96, EGM2008,
		   mean sea level, MSL,
		   height above ellipsoid, HAE,
		   vertical datum,
		   latitude and longitude,
		   online calculator,
		   WGS84 ellipsoid,
		   GeographicLib" />
  </head>
  <body>
    <h3>
      Online geoid calculations using the
      <a href="http://geographiclib.sourceforge.net/html/utilities.html#geoideval">
	GeoidEval</a> utility
    </h3>
    <form action="/cgi-bin/GeoidEval" method="get">
      <p>
        Position (ex. &laquo;<tt>16.78 -3.01</tt>&raquo;, &laquo;<tt>16d46'33"N 3d0.6'W</tt>&raquo;):<br>
        &nbsp;&nbsp;&nbsp;
        <input type=text name="input" size=30 value="$INPUTENC">
      </p>
      <p>
        Select action:<br>
        &nbsp;&nbsp;&nbsp;
        <input type="submit" name="option" value="Submit">
        <input type="submit" name="option" value="Reset">
      </p>
      <p>
        Geoid height:
<font size="4"><pre>
    lat lon = $POSITION1 `encodevalue "$POSITION2"`
    geoid heights (m)
	<a href="http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008">EGM2008</a> = $F`encodevalue "$HEIGHT2008"`$G
	<a href="http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html">EGM96</a>   = $F`encodevalue "$HEIGHT96"`$G
	<a href="http://earth-info.nga.mil/GandG/wgs84/gravitymod/wgs84_180/wgs84_180.html">EGM84</a>   = $F`encodevalue "$HEIGHT84"`$G</pre></font>
      </p>
    </form>
    <hr>
    <p>
      <a href="http://geographiclib.sourceforge.net/html/utilities.html#geoideval">
        GeoidEval</a>
      computes the height of the geoid above the WGS84 ellipsoid
      using interpolation in a grid of values for the earth
      gravity models,
      <a href="http://earth-info.nga.mil/GandG/wgs84/gravitymod/wgs84_180/wgs84_180.html">
        EGM84</a>, or
      <a href="http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html">
        EGM96</a>,
      <a href="http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008">
            EGM2008</a>.
      The RMS error in the interpolated height is about 1 mm.
      Give the position in terms of latitude and longitude, for example
      (these all refer to the position of Timbuktu):
      <pre>
        16.776 -3.009
        16d47' -3d1'
        W3d0'34" N16d46'33"</pre>
    </p>
    <p>
      <a href="http://geographiclib.sourceforge.net/html/utilities.html#geoideval">
        GeoidEval</a>,
      which is a simple wrapper of the
      <a href="http://geographiclib.sourceforge.net/html/classGeographicLib_1_1Geoid.html">
        GeographicLib::Geoid</a> class,
      is one of the utilities provided
      with <a href="http://geographiclib.sourceforge.net/">
        GeographicLib</a>.
      This web interface illustrates a subset of its capabilities.  If
      you wish to use GeoidEval directly,
      <a href="http://sourceforge.net/projects/geographiclib/files/distrib">
        download</a>
      and compile GeographicLib.  A description of the methods is given
      <a href="http://geographiclib.sourceforge.net/html/geoid.html">
        here</a>.
    </p>
    <p>
    </p>
    <hr>
    <address><a href="http://charles.karney.info/">Charles Karney</a>
      <a href="mailto:charles@karney.com">&lt;charles@karney.com&gt;</a>
      (2009-10-27)</address>
    <a href="http://sourceforge.net">
      <img src="http://sourceforge.net/sflogo.php?group_id=283628&amp;type=2" border="0" alt="SourceForge.net" />
    </a>
  </body>
</html>
EOF
