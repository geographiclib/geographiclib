#! /bin/sh
#
# tar.gz and zip distrib files copied to $DEVELSOURCE
# html documentation rsync'ed to  $DEVELSOURCE/doc/html/
#
# Windows version ready to build in
# $WINDOWSBUILD/GeographicLib-$VERSION/BUILD-vc10
# after compile + package installer is GeographicLib-$VERSION-win32.exe
#
# Built version ready to install in /usr/local in
# relc/GeographicLib-$VERSION/BUILD-system
#
# python update
# in gita/geographiclib/python ...
# python setup.py sdist --formats gztar,zip upload
#
# $Id$

VERSION=1.18
BRANCH=master
TEMP=/scratch/geographic-dist
GITSOURCE=file:///home/ckarney/afs/geographiclib
DEVELSOURCE=$HOME/geographiclib
WINDOWSBUILD=$HOME/afs/temp
test -d $TEMP || mkdir $TEMP
rm -rf $TEMP/*
mkdir $TEMP/gita
cd $TEMP/gita
git clone -b $BRANCH $GITSOURCE
cd geographiclib
sh autogen.sh
mkdir BUILD
cd BUILD
cmake ..
make dist
cp GeographicLib-$VERSION.{zip,tar.gz} $DEVELSOURCE
rsync -a --delete doc/html/ $DEVELSOURCE/doc/html/
mkdir $TEMP/rel{a,b,c,x,y}
tar xfpzC GeographicLib-$VERSION.tar.gz $TEMP/rela
tar xfpzC GeographicLib-$VERSION.tar.gz $TEMP/relb
tar xfpzC GeographicLib-$VERSION.tar.gz $TEMP/relc
tar xfpzC GeographicLib-$VERSION.tar.gz $TEMP/relx
rm -rf $WINDOWSBUILD/GeographicLib-$VERSION
tar xfpzC GeographicLib-$VERSION.tar.gz $WINDOWSBUILD
mkdir $WINDOWSBUILD/GeographicLib-$VERSION/BUILD-vc10
echo cmake -G \'Visual Studio 10\' -D ENABLE_MATLAB=ON -D CMAKE_INSTALL_PREFIX=C:/pkg-vc10/GeographicLib .. > $WINDOWSBUILD/GeographicLib-$VERSION/BUILD-vc10/config

cd $TEMP/rela/GeographicLib-$VERSION
make -j10
make PREFIX=$TEMP/insta install
cd $TEMP/insta
find . -type f | sort -u > ../files.a

cd $TEMP/relb/GeographicLib-$VERSION
./configure --prefix=$TEMP/instb
make -j10
make install
mv $TEMP/instb/share/doc/{geographiclib,GeographicLib}
cd $TEMP/instb
find . -type f | sort -u > ../files.b

cd $TEMP/relc/GeographicLib-$VERSION
mkdir BUILD
cd BUILD
cmake -D MAINTAINER=OFF -D CMAKE_INSTALL_PREFIX=$TEMP/instc ..
make -j10
make install
mkdir ../BUILD-maint
cd ../BUILD-maint
cmake -D MAINTAINER=ON -D CMAKE_INSTALL_PREFIX=$TEMP/instd ..
make -j10
make install
mkdir ../BUILD-matlab
cd ../BUILD-matlab
cmake -D MAINTAINER=OFF -D ENABLE_MATLAB=ON -D CMAKE_INSTALL_PREFIX=$TEMP/inste ..
make -j10
make -j10 matlab-all
make install
mkdir ../BUILD-system
cd ../BUILD-system
cmake -D MAINTAINER=OFF -D ENABLE_MATLAB=ON ..
make -j10
make -j10 matlab-all

cd $TEMP/instc
find . -type f | sort -u > ../files.c
cd $TEMP/instd
find . -type f | sort -u > ../files.d
cd $TEMP/inste
find . -type f | sort -u > ../files.e

mkdir $TEMP/gitb
cd $TEMP/gitb
git clone -b $BRANCH $GITSOURCE
cd geographiclib
sh autogen.sh
./configure
make dist-gzip
tar xfpzC geographiclib-$VERSION.tar.gz $TEMP/rely
mv $TEMP/rely/{geographiclib,GeographicLib}-$VERSION
cd $TEMP/rely
find . -type f | sort -u > ../files.y
cd $TEMP/relx
find . -type f | sort -u > ../files.x

cd $TEMP
cat > testprogram.cpp <<EOF
#include <iostream>
#include <iomanip>

#include <GeographicLib/Constants.hpp>
#include <GeographicLib/DMS.hpp>
#include <GeographicLib/LambertConformalConic.hpp>

int main() {
  using namespace GeographicLib;
  double 
    // These are the constants for Pennsylvania South, EPSG:3364
    // http://www.spatialreference.org/ref/epsg/3364/
    a = Constants::WGS84_a(),   // major radius
    r = 298.257222101,          // inverse flattening (GRS80)
    lat1 = DMS::Decode(40,58),  // standard parallel 1
    lat2 = DMS::Decode(39,56),  // standard parallel 2
    k1 = 1,                     // scale on std parallels
    lat0 =  DMS::Decode(39,20), // latitude of origin
    lon0 = -DMS::Decode(77,45), // longitude of origin
    fe = 600000,                // false easting
    fn = 0;                     // false northing
  LambertConformalConic PASouth(a, r, lat1, lat2, k1);
  double x0, y0;
  PASouth.Forward(lon0, lat0, lon0, x0, y0); // Transform origin point
  x0 -= fe; y0 -= fn;           // Combine result with false origin

  double lat = 39.95, lon = -75.17;    // Philadelphia
  double x, y;
  PASouth.Forward(lon0, lat, lon, x, y);
  x -= x0; y -= y0;             // Philadelphia in PA South coordinates

  std::cout << std::fixed << std::setprecision(3)
            << x << " " << y << "\n";
  return 0;
}
EOF
for i in a b c d e; do
    cp testprogram.cpp testprogram$i.cpp
    g++ -c -g -O3 -I$TEMP/inst$i/include testprogram$i.cpp
    g++ -g -o testprogram$i testprogram$i.o -Wl,-rpath=$TEMP/inst$i/lib -L$TEMP/inst$i/lib -lGeographic
    ./testprogram$i
done
