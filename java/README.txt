This is a Java implementation of the geodesic algorithms described in

  C. F. F. Karney,
  Algorithms for geodesics,
  J. Geodesy 87, 43-55 (2013);
  http://dx.doi.org/10.1007/s00190-012-0578-z
  Addenda: http://geographiclib.sf.net/geod-addenda.html

For documentation, see

  http://geographiclib.sourceforge.net/html/java/

The code in this directory is entirely self-contained.  In particular,
it does not depend on the C++ classes.  You can build the example
programs using, for example,

  cd examples
  javac -cp .:../src/main/java Inverse.java
  echo -30 0 29.5 179.5 | java -cp .:../src/main/java Inverse

Alternatively you can create a jar file for the libary and use that

  mvn package
  cd examples
  javac -cp .:../targets/GeographicLib-1.31.jar Inverse.java
  echo -30 0 29.5 179.5 |
    java -cp .:../targets/GeographicLib-1.31.jar Inverse
