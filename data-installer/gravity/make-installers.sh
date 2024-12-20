#! /bin/sh
for m in egm84 egm96 egm2008 wgs84; do
  tar cfCjv ../../data-distrib/gravity-distrib/$m.tar.bz2 .. gravity/$m.egm{,.cof}
  (
      cd ..
      rm -f ../data-distrib/gravity-distrib/$m.zip
      zip ../data-distrib/gravity-distrib/$m.zip gravity/$m.egm{,.cof}
  )
done
