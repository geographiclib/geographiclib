#! /bin/sh
#
# Planimeter.cgi
# cgi script for measuring the area of geodesic polygons
#
# Copyright (c) Charles Karney (2011-2019) <charles@karney.com> and
# licensed under the MIT/X11 License.  For more information, see
# https://geographiclib.sourceforge.io/

. ./utils.sh
OPTION=`lookupkey "$QUERY_STRING" option`
if test "$OPTION" = Reset; then
    INPUT=
else
    INPUT=`lookupcheckkey "$QUERY_STRING" input`
    NORM=`lookupkey "$QUERY_STRING" norm`
    TYPE=`lookupkey "$QUERY_STRING" type`
    RHUMB=`lookupkey "$QUERY_STRING" rhumb`
fi
env > /tmp/env

INPUTENC=`encodevalue "$INPUT"`
if test "$TYPE" = "polyline"; then
  LINEFLAG=-l
else
  LINEFLAG=
fi
if test "$RHUMB" = "rhumb"; then
  RHUMBFLAG=-R
else
  RHUMBFLAG=
fi
EXECDIR=../bin
COMMAND="Planimeter"
VERSION=`$EXECDIR/$COMMAND --version | cut -f4 -d" "`
STATUS=
NUM=
LEN=
AREA=
if test "$INPUT"; then
    STATUS=`echo "$INPUT" | head -1 | $EXECDIR/GeoConvert`
    if test $? -eq 0; then
        STATUS=OK
        OUTPUT=`echo "$INPUT" | $EXECDIR/$COMMAND $LINEFLAG $RHUMBFLAG |
            head -1`
        NUM="`echo $OUTPUT | cut -f1 -d' '`"
        LEN="`echo $OUTPUT | cut -f2 -d' '`"
        AREA="`echo $OUTPUT | cut -f3 -d' '`"
        if test "$NORM"; then
            TRANSFORMEDINPUT=`echo "$INPUT" | $EXECDIR/GeoConvert -p 20 |
            sed '/^ERROR/,$d'`
            INPUTENC=`encodevalue "$TRANSFORMEDINPUT"`
        fi
    fi
    # echo `date +"%F %T"` "$COMMAND: $INPUT" >> ../persistent/utilities.log
fi

echo Content-type: text/html
echo
cat <<EOF
<html>
  <head>
    <title>
      Online geodesic planimeter
    </title>
    <meta name="description" content="Online geodesic planimeter" />
    <meta name="author" content="Charles F. F. Karney" />
    <meta name="keywords"
          content="geodesics,
                   geodesic distance,
                   geodesic area,
                   geographic distance,
                   geographic area,
                   geodesic polygon,
                   shortest path,
                   spheroidal triangle,
                   latitude and longitude,
                   rhumb lines,
                   online calculator,
                   WGS84 ellipsoid,
                   GeographicLib" />
  </head>
  <body>
    <h3>
      Online geodesic polygon calculations using the
      <a href="https://geographiclib.sourceforge.io/html/Planimeter.1.html">
         Planimeter</a> utility
    </h3>
    <form action="/cgi-bin/Planimeter" method="get">
      <p>
        <table>
          <tr>
            <td>
              Closed/open:
            </td>
            <td>&nbsp;<label for="c">
                <input type="radio" name="type" value="polygon" id="c"
                       `test "$LINEFLAG" || echo CHECKED`>
                Polygon</label>
            </td>
            <td>&nbsp;<label for="o">
                <input type="radio" name="type" value="polyline" id="o"
                       `test "$LINEFLAG" && echo CHECKED`>
                Polyline</label>
            </td>
          </tr>
          <tr>
            <td>
              Edge type:
            </td>
            <td>&nbsp;<label for="g">
                <input type="radio" name="rhumb" value="geodesic" id="g"
                       `test "$RHUMBFLAG" || echo CHECKED`>
                Geodesic</label>
            </td>
            <td>&nbsp;<label for="r">
                <input type="radio" name="rhumb" value="rhumb" id="r"
                       `test "$RHUMBFLAG" && echo CHECKED`>
                Rhumb line</label>
            </td>
          </tr>
        </table>
      </p>
      <p>
        Enter the vertices as latitude longitude pairs, one per line:
        <br>
        &nbsp;&nbsp;&nbsp;
        <textarea cols=45 rows=15 name="input">$INPUTENC</textarea>
      </p>
      <p>
        <input type="checkbox" name="norm" value="decdegrees"
               `test "$NORM" && echo CHECKED` />
        Convert vertices to decimal degrees
      </p>
      <p>
        Select action:<br>
        &nbsp;&nbsp;&nbsp;
        <input type="submit" name="option" value="Submit">
        <input type="submit" name="option" value="Reset">
      </p>
      <p>
        Results:<br>
        <font size="4"><pre>
    STATUS             = $STATUS
    number of vertices = $NUM
    Perimeter (m)      = $LEN
    area (m^2)         = $AREA</pre></font>
      </p>
    </form>
    <hr>
    <p>
      In polygon mode,
      <a href="https://geographiclib.sourceforge.io/html/Planimeter.1.html">
        Planimeter (version $VERSION)</a>
      calculates the perimeter and area of a polygon whose edges are
      either geodesics or rhumb lines on the WGS84 ellipsoid.
      Counter-clockwise traversal of a polygon results in a positive
      area.  Arbitrarily complex polygons are allowed.  In the case of
      self-intersecting polygons the area is accumulated
      "algebraically", i.e., the areas of the 2 loops in a figure-8
      polygon will partially cancel.  There is no need to close the
      polygon.  Polygons may include one or both poles.  In polyline
      mode,
      <a href="https://geographiclib.sourceforge.io/html/Planimeter.1.html">
        Planimeter</a>
      calculates the length of the geodesic path joining the points.
    </p>
    <p>
      Give the vertices in terms of latitude and longitude, for example
      (these all refer to the position of Timbuktu):
      <pre>
        16.776 -3.009
        16d47' -3d1'
        W3&deg;0'34" N16&deg;46'33"
        3:0:34W 16:46:33N</pre>
      The coordinates can also be given in UTM, UPS, or MGRS coordinates (see
      the documentation on the
      <a href="https://geographiclib.sourceforge.io/html/GeoConvert.1.html">
        GeoConvert</a>
      utility).  A blank line or a coordinate which cannot be understood
      causes the reading of vertices to be stopped.
    </p>
    <p>
      The result for the perimeter is accurate to about 15&nbsp;nm per
      vertex.  The result for the area is accurate to about
      0.1&nbsp;m<sup>2</sup> per vertex.
    </p>
    <p>
      <a href="https://geographiclib.sourceforge.io/html/Planimeter.1.html">
        Planimeter</a>,
      which is a simple wrapper of the
      <a href="https://geographiclib.sourceforge.io/html/classGeographicLib_1_1PolygonAreaT.html">
        GeographicLib::PolygonAreaT</a> class,
      is one of the utilities provided
      with <a href="https://geographiclib.sourceforge.io/">
        GeographicLib</a>.
      Geodesic areas can also be computed using JavaScript; see the
      <a href="../scripts/geod-calc.html">JavaScript geodesic
        calculator</a>.
      If you wish to use Planimeter directly,
      <a href="https://sourceforge.net/projects/geographiclib/files/distrib">
        download</a>
      and compile GeographicLib.  The algorithms are described
      in C. F. F. Karney,
      <a href="https://doi.org/10.1007/s00190-012-0578-z"><i>Algorithms for
      geodesics</i></a>,
      J. Geodesy <b>87</b>, 43&ndash;55 (2013); DOI:
      <a href="https://doi.org/10.1007/s00190-012-0578-z">
        10.1007/s00190-012-0578-z</a>;
      addenda:
      <a href="https://geographiclib.sourceforge.io/geod-addenda.html">
        geod-addenda.html</a>.
    </p>
    <hr>
    <address>Charles Karney
      <a href="mailto:charles@karney.com">&lt;charles@karney.com&gt;</a>
      (2015-02-24)</address>
    <a href="https://geographiclib.sourceforge.io">
      GeographicLib home
    </a>
  </body>
</html>
EOF
