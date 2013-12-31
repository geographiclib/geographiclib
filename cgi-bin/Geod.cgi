#! /bin/sh
#
# Geod.cgi (redirect to GeodSolve.cgi)
#
# Copyright (c) Charles Karney (2013) <charles@karney.com> and licensed
# under the MIT/X11 License.  For more information, see
# http://geographiclib.sourceforge.net/

cat <<EOF
Content-type: text/html

<html>
  <head>
    <title>Online geodesic calculator</title>
    <meta HTTP-EQUIV="Refresh"
	  CONTENT="5; URL=http://geographiclib.sourceforge.net/cgi-bin/GeodSolve">
  </head>
  <body topmargin=10 leftmargin=10>
    <h3>
      <blockquote>
	<em>
	  The Geod calculator has been renamed GeodSolve which is available at
	  <center>
	    <a href="http://geographiclib.sourceforge.net/cgi-bin/GeodSolve">
	      http://geographiclib.sourceforge.net/cgi-bin/GeodSolve</a>.
	  </center>
	  <br>
	  You will be redirected there.  Click on the link to go there
	  directly.
	</em>
      </blockquote>
    </h3>
  </body>
</html>
EOF
