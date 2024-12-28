#! /bin/sh
#for m in emm2010 wmm2010 igrf11 wmm2015 wmm2015v2 wmm2020 igrf13 wmm2025 wmmhr2025 igrf14; do
for m in igrf14; do
  tar cfCjv ../../data-distrib/magnetic-distrib/$m.tar.bz2 .. magnetic/$m.wmm{,.cof}
  (
      cd ..
      rm -f ../data-distrib/magnetic-distrib/$m.zip
      zip ../data-distrib/magnetic-distrib/$m.zip magnetic/$m.wmm{,.cof}
  )
done
