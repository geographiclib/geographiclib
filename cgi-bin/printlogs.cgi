#!/bin/sh
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
