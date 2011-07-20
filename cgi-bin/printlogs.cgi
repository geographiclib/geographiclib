#! /bin/sh
#
# printlogs.cgi
# cgi script for printing error logs
#
# Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed
# under the MIT/X11 License.  For more information, see
# http://geographiclib.sourceforge.net/
#
# $Id$
echo Content-type: text/html
echo
cat <<EOF
<html>
  <header>
    <title>
      GeographicLib web utilities log
    </title>
  </header>
  <body>
    <pre>
EOF
cat ../persistent/utilities.log
cat <<EOF
    </pre>
  </body>
</html>
