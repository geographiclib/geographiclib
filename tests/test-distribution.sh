#! /bin/sh
#
# tar.gz and zip distrib files copied to $DEVELSOURCE
# html documentation rsync'ed to  $DEVELSOURCE/doc/html/
#
# Windows version ready to build in
# $WINDOWSBUILD/GeographicLib-$VERSION/BUILD-vc10{,-x64}
# after ./build installer is copied to
# $DEVELSOURCE/GeographicLib-$VERSION-win{32,64}.exe
#
# Built version ready to install in /usr/local in
# relc/GeographicLib-$VERSION/BUILD-system
#
# python update
# in gita/geographiclib/python ...
# python setup.py sdist --formats gztar,zip upload
#
# gita - check out from git, create package with cmake
# gita - check out from git, create package with autoconf
# rela - release package, build with make
# relb - release package, build with autoconf
# relc - release package, build with cmake
# relx - cmake release package inventory
# rely - autoconf release package inventory
# insta - installed files, make
# instb - installed files, autoconf
# instc - installed files, cmake
# inste - installed files, cmake matlab=on
# instf - installed files, autoconf direct from git repository

set -e

VERSION=1.23
BRANCH=devel
TEMP=/scratch/geographic-dist
DEVELSOURCE=/u/geographiclib
GITSOURCE=file://$DEVELSOURCE
WINDOWSBUILD=/u/temp

test -d $TEMP || mkdir $TEMP
rm -rf $TEMP/*
mkdir $TEMP/gita # Package creation via cmake
mkdir $TEMP/gitb # Package creation via autoconf
(cd $TEMP/gita; git clone -b $BRANCH $GITSOURCE)
(cd $TEMP/gitb; git clone -b $BRANCH $GITSOURCE)
cd $TEMP/gita/geographiclib
sh autogen.sh
mkdir BUILD
cd BUILD
cmake ..
make dist
cp GeographicLib-$VERSION.{zip,tar.gz} $DEVELSOURCE
rsync -a --delete doc/html/ $DEVELSOURCE/doc/html/

mkdir $TEMP/rel{a,b,c,x,y}
tar xfpzC GeographicLib-$VERSION.tar.gz $TEMP/rela # Version of make build
tar xfpzC GeographicLib-$VERSION.tar.gz $TEMP/relb # Version for autoconf
tar xfpzC GeographicLib-$VERSION.tar.gz $TEMP/relc # Version for cmake
tar xfpzC GeographicLib-$VERSION.tar.gz $TEMP/relx
rm -rf $WINDOWSBUILD/GeographicLib-$VERSION
unzip -qq -d $WINDOWSBUILD GeographicLib-$VERSION.zip 
mkdir $WINDOWSBUILD/GeographicLib-$VERSION/BUILD-vc10
(
    echo cmake -G \"Visual Studio 10\" -D PACKAGE_PATH=u:/pkg-vc10 -D GEOGRAPHICLIB_EXAMPLES=ON ..
    echo cmake --build . --config Release --target ALL_BUILD
    echo cmake --build . --config Release --target RUN_TESTS
    echo cmake --build . --config Release --target INSTALL
    echo cmake --build . --config Release --target PACKAGE
    echo cp GeographicLib-$VERSION-win32.exe $DEVELSOURCE/
) > $WINDOWSBUILD/GeographicLib-$VERSION/BUILD-vc10/build
chmod +x $WINDOWSBUILD/GeographicLib-$VERSION/BUILD-vc10/build
mkdir $WINDOWSBUILD/GeographicLib-$VERSION/BUILD-vc10-shared
(
    echo cmake -G \"Visual Studio 10\" -D PACKAGE_PATH=u:/pkg-vc10-shared -D GEOGRAPHICLIB_EXAMPLES=ON -D GEOGRAPHIC_SHARED_LIB=ON ..
    echo cmake --build . --config Release --target ALL_BUILD
    echo cmake --build . --config Release --target RUN_TESTS
    echo cmake --build . --config Release --target INSTALL
) > $WINDOWSBUILD/GeographicLib-$VERSION/BUILD-vc10-shared/build
chmod +x $WINDOWSBUILD/GeographicLib-$VERSION/BUILD-vc10-shared/build
mkdir $WINDOWSBUILD/GeographicLib-$VERSION/BUILD-vc10-x64
(
    echo cmake -G \"Visual Studio 10 Win64\" -D PACKAGE_PATH=u:/pkg-vc10-x64 -D GEOGRAPHICLIB_EXAMPLES=ON -D MATLAB_COMPILER=mex ..
    echo cmake --build . --config Release --target ALL_BUILD
    echo cmake --build . --config Release --target matlab-all
    echo cmake --build . --config Release --target RUN_TESTS
    echo cmake --build . --config Release --target INSTALL
    echo cmake --build . --config Release --target PACKAGE
    echo cp GeographicLib-$VERSION-win64.exe $DEVELSOURCE/
) > $WINDOWSBUILD/GeographicLib-$VERSION/BUILD-vc10-x64/build
chmod +x $WINDOWSBUILD/GeographicLib-$VERSION/BUILD-vc10-x64/build

cd $TEMP/rela/GeographicLib-$VERSION
make -j10
make PREFIX=$TEMP/insta install
cd $TEMP/insta
find . -type f | sort -u > ../files.a

cd $TEMP/relb/GeographicLib-$VERSION
mkdir BUILD-config
cd BUILD-config
../configure --prefix=$TEMP/instb
make -j10
make install
mv $TEMP/instb/share/doc/{geographiclib,GeographicLib}
cd $TEMP/instb
find . -type f | sort -u > ../files.b

cd $TEMP/relc/GeographicLib-$VERSION
mkdir BUILD
cd BUILD
cmake -D CMAKE_INSTALL_PREFIX=$TEMP/instc ..
make -j10
make install
mkdir ../BUILD-matlab
cd ../BUILD-matlab
cmake -D MATLAB_COMPILER=mkoctfile -D CMAKE_INSTALL_PREFIX=$TEMP/inste ..
make -j10
make -j10 matlab-all
make install
mkdir ../BUILD-system
cd ../BUILD-system
cmake -D MATLAB_COMPILER=mkoctfile ..
make -j10
make -j10 matlab-all

cd $TEMP/instc
find . -type f | sort -u > ../files.c
cd $TEMP/inste
find . -type f | sort -u > ../files.e

cd $TEMP/gitb/geographiclib
sh autogen.sh
mkdir BUILD-config
cd BUILD-config
../configure --prefix=$TEMP/instf
make dist-gzip
make install
tar xfpzC geographiclib-$VERSION.tar.gz $TEMP/rely
mv $TEMP/rely/{geographiclib,GeographicLib}-$VERSION
cd $TEMP/rely
find . -type f | sort -u > ../files.y
cd $TEMP/relx
find . -type f | sort -u > ../files.x

mv $TEMP/instf/share/doc/{geographiclib,GeographicLib}
cd $TEMP/instf
find . -type f | sort -u > ../files.f


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
for i in a b c d e f; do
    cp testprogram.cpp testprogram$i.cpp
    g++ -c -g -O3 -I$TEMP/inst$i/include testprogram$i.cpp
    g++ -g -o testprogram$i testprogram$i.o -Wl,-rpath=$TEMP/inst$i/lib -L$TEMP/inst$i/lib -lGeographic
    ./testprogram$i
done
