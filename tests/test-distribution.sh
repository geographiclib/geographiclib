#! /bin/sh
#
# tar.gz and zip distrib files copied to $DEVELSOURCE
# html documentation rsync'ed to $WEBDIST/htdocs/$VERSION-pre/
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
# [or: python setup.py sdist --formats gztar,zip register upload]
#
# gita - check out from git, create package with cmake
# gitb - check out from git, create package with autoconf
# gitr - new release branch
# rela - release package, build with make
# relb - release package, build with autoconf
# relc - release package, build with cmake
# relx - cmake release package inventory
# rely - autoconf release package inventory
# insta - installed files, make
# instb - installed files, autoconf
# instc - installed files, cmake
# instf - installed files, autoconf direct from git repository

set -e
umask 0022

# The following files contain version information:
#   pom.xml
#   CMakeLists.txt (PROJECT_VERSION_* LIBVERSION_*)
#   NEWS
#   configure.ac (AC_INIT, GEOGRAPHICLIB_VERSION_* LT_*)
#   tests/test-distribution.sh
#   doc/GeographicLib.dox.in (3 places)
#   doc/NETGeographicLib.dox (a few places)

# Need updating if underlying library changes

# python
#   python/setup.py
#   python/geographiclib/__init__.py
#   python/doc/index.rst (date + update change log)
#   python/README.rst
# use: cd python; pychecker geographiclib/*.py

# MATLAB
#   matlab/geographiclib/Contents.m version (multiple places) + date
#   matlab/geographiclib-blurb.txt version (multiple places) + date
#   update version number "%2F15%2F" in documentation link in index.html,
#     GeographicLib.dox.in, geodesic-{c,for}.dox,
#     java/src/main/java/net/sf/geographiclib/package-info.java,
#     js/GeographicLib.md, python/doc/index.rst
#   mathworks has switched to an uglier URL.  Only update if there are
#   changes.
# use MATLAB's analyze code

# C
#   legacy/C/geodesic.h comment + GEODESIC_VERSION_*
#   doc/geodesic-c.dox (date + update change log)
# PROJ integration
#   geodesic.[hc], geodtest.c (renamed to geodtest.cpp)
#   plus man/man1/geod.1 man/man3/geodesic.3

# Fortran
#   legacy/Fortran/geodesic.for comment + geover
#   doc/geodesic-for.dox (date + update change log)

# Java
#   java/pom.xml java/*/pom.xml
#   java/src/main/java/net/sf/geographiclib/package-info.java (date +
#   update change log)
#   (remember to remove SNAPSHOT from version number of lib in test programs)
# probably should deploy a SNAPSHOT version of the lib to grease the
# wheels

# maxima
#   maxima/geodesic.mac

# JavaScript
#   js/src/Math.js
#   js/package.json
#   js/README.md
#   js/GeographicLib.md (date + update change log)
# use: cd js; jshint src

DATE=`date +%F`
VERSION=1.50.1
BRANCH=devel
TEMP=/scratch/geographiclib-dist
if test `hostname` = petrel.petrel.org; then
    DEVELSOURCE=$HOME/geographiclib
    WINDEVELSOURCE=/w/geographiclib
    WINDOWSBUILD=/var/tmp
else
    DEVELSOURCE=/u/geographiclib
    WINDEVELSOURCE=/w/geographiclib
    WINDOWSBUILD=/u/temp
fi
WINDOWSBUILDWIN=w:/temp
GITSOURCE=file://$DEVELSOURCE
WEBDIST=/home/ckarney/web/geographiclib-web
NUMCPUS=4
HAVEINTEL=

test -d $TEMP || mkdir $TEMP
rm -rf $TEMP/*
mkdir $TEMP/gita # Package creation via cmake
mkdir $TEMP/gitb # Package creation via autoconf
mkdir $TEMP/gitr # For release branch
(cd $TEMP/gitr; git clone -b $BRANCH $GITSOURCE)
(cd $TEMP/gita; git clone -b $BRANCH file://$TEMP/gitr/geographiclib)
(cd $TEMP/gitb; git clone -b $BRANCH file://$TEMP/gitr/geographiclib)
cd $TEMP/gita/geographiclib
sh autogen.sh
mkdir BUILD
cd BUILD
cmake -D GEOGRAPHICLIB_LIB_TYPE=BOTH -D GEOGRAPHICLIB_DOCUMENTATION=ON ..
make dist
cp GeographicLib-$VERSION.{zip,tar.gz} $DEVELSOURCE
make doc
(
    cd ../java
    mvn -q package -P release
    rsync -a target/apidocs/ ../BUILD/doc/html/java/
)
(
    cd ../python
    python2 -m unittest -v geographiclib.test.test_geodesic
    python3 -m unittest -v geographiclib.test.test_geodesic
)
(
    cd ../matlab/geographiclib
    octave --no-gui --no-window-system --eval geographiclib_test
)
(
   cd js/geographiclib
   npm test
)
rsync -a --delete doc/html/ $WEBDIST/htdocs/$VERSION-pre/
mkdir -p $TEMP/js
cp -p js/*.js js/*.html $TEMP/js/
JS_VERSION=`grep Version: $TEMP/js/geographiclib.js | cut -f2 -d: | tr -d ' '`
mv $TEMP/js/geographiclib.js $TEMP/js/geographiclib-$JS_VERSION.js
ln -s geographiclib-$JS_VERSION.js $TEMP/js/geographiclib.js
mv $TEMP/js/geographiclib.min.js $TEMP/js/geographiclib-$JS_VERSION.min.js
ln -s geographiclib-$JS_VERSION.min.js $TEMP/js/geographiclib.min.js
rsync -a --delete $TEMP/js/ $WEBDIST/htdocs/scripts/test/

mkdir $TEMP/rel{a,b,c,x,y}
tar xfpzC GeographicLib-$VERSION.tar.gz $TEMP/rela # Version of make build
tar xfpzC GeographicLib-$VERSION.tar.gz $TEMP/relb # Version for autoconf
tar xfpzC GeographicLib-$VERSION.tar.gz $TEMP/relc # Version for cmake + mvn
tar xfpzC GeographicLib-$VERSION.tar.gz $TEMP/relx
rm -rf $WINDOWSBUILD/GeographicLib-$VERSION
unzip -qq -d $WINDOWSBUILD GeographicLib-$VERSION.zip

cat > $WINDOWSBUILD/GeographicLib-$VERSION/mvn-build <<'EOF'
#! /bin/sh -exv
unset GEOGRAPHICLIB_DATA
# for v in 2019 2017 2015 2013 2012 2010; do
for v in 2019 2017 2015; do
  for a in 64 32; do
    echo ========== maven $v-$a ==========
    rm -rf c:/scratch/geog-mvn-$v-$a
    mvn -Dcmake.compiler=vc$v -Dcmake.arch=$a \
      -Dcmake.project.bin.directory=c:/scratch/geog-mvn-$v-$a install
  done
done
EOF
chmod +x $WINDOWSBUILD/GeographicLib-$VERSION/mvn-build
cp $TEMP/gita/geographiclib/pom.xml $WINDOWSBUILD/GeographicLib-$VERSION/

# for ver in 10 11 12 14 15 16; do
for ver in 14 15 16; do
    for arch in win32 x64; do
	pkg=vc$ver-$arch
	gen="Visual Studio $ver"
	installer=
	# N.B. update CPACK_NSIS_INSTALL_ROOT in CMakeLists.txt and
	# update documentation examples if VS version for binary
	# installer changes.
	test "$ver" = 14 && installer=y
	mkdir $WINDOWSBUILD/GeographicLib-$VERSION/BUILD-$pkg
	(
	    echo "#! /bin/sh -exv"
	    echo echo ========== cmake $pkg ==========
	    echo 'b=geog-`pwd | sed s%.*/%%`'
	    echo rm -rf c:/scratch/\$b w:/pkg-$pkg/GeographicLib-$VERSION/\*
	    echo 'unset GEOGRAPHICLIB_DATA'
	    echo 'mkdir -p c:/scratch/$b'
	    echo 'cd c:/scratch/$b'
	    echo cmake -G \"$gen\" -A $arch -D GEOGRAPHICLIB_LIB_TYPE=BOTH -D CMAKE_INSTALL_PREFIX=w:/pkg-$pkg/GeographicLib-$VERSION -D PACKAGE_DEBUG_LIBS=ON -D BUILD_NETGEOGRAPHICLIB=ON -D CONVERT_WARNINGS_TO_ERRORS=ON $WINDOWSBUILDWIN/GeographicLib-$VERSION
	    echo cmake --build . --config Debug   --target ALL_BUILD
	    echo cmake --build . --config Debug   --target RUN_TESTS
	    echo cmake --build . --config Debug   --target INSTALL
	    echo cmake --build . --config Release --target ALL_BUILD
	    echo cmake --build . --config Release --target exampleprograms
	    echo cmake --build . --config Release --target netexamples
	    echo cmake --build . --config Release --target RUN_TESTS
	    echo cmake --build . --config Release --target INSTALL
	    echo cmake --build . --config Release --target PACKAGE
	    test "$installer" &&
		echo cp GeographicLib-$VERSION-*.exe $WINDEVELSOURCE/ || true
	    echo 'b=geogc-`pwd | sed s%.*/%%`'
	    echo rm -rf c:/scratch/\$b
	    echo 'mkdir -p c:/scratch/$b'
	    echo 'cd c:/scratch/$b'
	    echo cmake -G \"$gen\" -A $arch -D CONVERT_WARNINGS_TO_ERRORS=ON $WINDOWSBUILDWIN/GeographicLib-$VERSION/legacy/C
	    echo cmake --build . --config Debug   --target ALL_BUILD
	    echo cmake --build . --config Debug   --target RUN_TESTS
	    echo cmake --build . --config Release --target ALL_BUILD
	    echo cmake --build . --config Release --target RUN_TESTS
	) > $WINDOWSBUILD/GeographicLib-$VERSION/BUILD-$pkg/build
	chmod +x $WINDOWSBUILD/GeographicLib-$VERSION/BUILD-$pkg/build
    done
done
cat > $WINDOWSBUILD/GeographicLib-$VERSION/test-all <<'EOF'
#! /bin/sh
(
    # Queue vs2015 build first for the binary installers
    for d in BUILD-vc14* BUILD-vc*; do
	test -f $d/build.done && continue
	(cd $d; ./build; touch build.done)
    done
    ./mvn-build
) >& build.log
EOF
chmod +x $WINDOWSBUILD/GeographicLib-$VERSION/test-all

cd $TEMP/gitr/geographiclib
git checkout release
git config user.email karney@users.sourceforge.net
find . -type f | grep -v '/\.git' | xargs rm
tar xfpz $DEVELSOURCE/GeographicLib-$VERSION.tar.gz
(
    cd GeographicLib-$VERSION
    find . -type f | while read f; do
	dest=../`dirname $f`
	test -d $dest || mkdir -p $dest
	mv $f $dest/
    done
)
rm -rf GeographicLib-$VERSION
find * -type d -empty | xargs -r rmdir
find * -type d -empty | xargs -r rmdir

cd $TEMP/rela/GeographicLib-$VERSION
make -j$NUMCPUS
make PREFIX=$TEMP/insta install
cd $TEMP/insta
find . -type f | sort -u > ../files.a

cd $TEMP/relb/GeographicLib-$VERSION
mkdir BUILD-config
cd BUILD-config
../configure --prefix=$TEMP/instb
make -j$NUMCPUS
make install
cd ..

if test "$HAVEINTEL"; then
    mkdir BUILD-config-intel
    cd BUILD-config-intel
    env FC=ifort CC=icc CXX=icpc ../configure
    make -j$NUMCPUS
    cd ..
fi

mv $TEMP/instb/share/doc/{geographiclib,GeographicLib}
cd $TEMP/instb
find . -type f | sort -u > ../files.b

cd $TEMP/relc/GeographicLib-$VERSION
mkdir BUILD
cd BUILD
cmake -D GEOGRAPHICLIB_LIB_TYPE=BOTH -D GEOGRAPHICLIB_DOCUMENTATION=ON -D USE_BOOST_FOR_EXAMPLES=ON -D CONVERT_WARNINGS_TO_ERRORS=ON -D CMAKE_INSTALL_PREFIX=$TEMP/instc ..
make -j$NUMCPUS all
make test
make -j$NUMCPUS exampleprograms
make install
(
    cd $TEMP/instc/lib/node_modules/geographiclib
    mocha
)

mkdir ../BUILD-system
cd ../BUILD-system
cmake -D GEOGRAPHICLIB_LIB_TYPE=BOTH -D CONVERT_WARNINGS_TO_ERRORS=ON ..
make -j$NUMCPUS all
make test
cd ..

if test "$HAVEINTEL"; then
    mkdir ../BUILD-intel
    cd ../BUILD-intel
    env FC=ifort CC=icc CXX=icpc cmake -D GEOGRAPHICLIB_LIB_TYPE=BOTH -D CONVERT_WARNINGS_TO_ERRORS=ON ..
    make -j$NUMCPUS all
    make test
    make -j$NUMCPUS exampleprograms
    cd ..
fi

# mvn -Dcmake.project.bin.directory=$TEMP/mvn install

cd $TEMP/gita/geographiclib/tests/sandbox
mkdir BUILD
cd BUILD
cmake -D CMAKE_PREFIX_PATH=$TEMP/instc ..
make

cd $TEMP/instc/share/matlab/geographiclib
mkdir $TEMP/matlab
cp -pr $TEMP/instc/share/matlab/geographiclib $TEMP/matlab
cd $TEMP/matlab/geographiclib
rm -f $DEVELSOURCE/geographiclib_toolbox_$VERSION.zip
zip $DEVELSOURCE/geographiclib_toolbox_$VERSION.zip *.m private/*.m
cd $TEMP/matlab
cp -p $TEMP/gita/geographiclib/geodesic.png .
cp -p $TEMP/gita/geographiclib/matlab/geographiclib-blurb.txt .
VERSION=$VERSION DATE=$DATE ROOT=$TEMP/matlab \
       sh $DEVELSOURCE/tests/matlab-toolbox-config.sh

cp -pr $TEMP/relc/GeographicLib-$VERSION/legacy $TEMP/
for l in C Fortran; do
    (
	mkdir $TEMP/legacy/$l/BUILD
	cd $TEMP/legacy/$l/BUILD
	cmake -D CONVERT_WARNINGS_TO_ERRORS=ON ..
	make -j$NUMCPUS all
	make test
	test $l = Fortran && continue
	if test "$HAVEINTEL"; then
	    mkdir $TEMP/legacy/$l/BUILD-intel
	    cd $TEMP/legacy/$l/BUILD-intel
	    env FC=ifort CC=icc CXX=icpc \
		cmake -D CONVERT_WARNINGS_TO_ERRORS=ON ..
	    make -j$NUMCPUS all
	    make test
	fi
    )
done

cd $TEMP/gita/geographiclib
(
    cd BUILD
    make -j$NUMCPUS testprograms
)
cp $DEVELSOURCE/include/mpreal.h include/
for p in 1 3 4 5; do
    mkdir BUILD-$p
    (
	cd BUILD-$p
	cmake -D USE_BOOST_FOR_EXAMPLES=ON -D GEOGRAPHICLIB_PRECISION=$p ..
	make -j$NUMCPUS all
	if test $p -ne 1; then
	    make test
	fi
	make -j$NUMCPUS testprograms
    )
done

cd $TEMP/instc
find . -type f | sort -u > ../files.c

cd $TEMP/gitb/geographiclib
./autogen.sh
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
    // https://www.spatialreference.org/ref/epsg/3364/
    a = Constants::WGS84_a(),   // major radius
    f = 1/298.257222101,        // inverse flattening (GRS80)
    lat1 = DMS::Decode(40,58),  // standard parallel 1
    lat2 = DMS::Decode(39,56),  // standard parallel 2
    k1 = 1,                     // scale on std parallels
    lat0 =  DMS::Decode(39,20), // latitude of origin
    lon0 = -DMS::Decode(77,45), // longitude of origin
    fe = 600000,                // false easting
    fn = 0;                     // false northing
  LambertConformalConic PASouth(a, f, lat1, lat2, k1);
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
for i in a b c f; do
    cp testprogram.cpp testprogram$i.cpp
    g++ -c -g -O3 -I$TEMP/inst$i/include testprogram$i.cpp
    g++ -g -o testprogram$i testprogram$i.o -Wl,-rpath=$TEMP/inst$i/lib -L$TEMP/inst$i/lib -lGeographic
    ./testprogram$i
done

libversion=`find $TEMP/instc/lib -type f \
-name 'libGeographic.so.*' -printf "%f" |
sed 's/libGeographic\.so\.//'`
test -f $TEMP/instb/lib/libGeographic.so.$libversion ||
echo autoconf/cmake library so mismatch

CONFIG_FILE=$TEMP/gitr/geographiclib/configure
CONFIG_MAJOR=`grep ^GEOGRAPHICLIB_VERSION_MAJOR= $CONFIG_FILE | cut -f2 -d=`
CONFIG_MINOR=`grep ^GEOGRAPHICLIB_VERSION_MINOR= $CONFIG_FILE | cut -f2 -d=`
CONFIG_PATCH=`grep ^GEOGRAPHICLIB_VERSION_PATCH= $CONFIG_FILE | cut -f2 -d=`
CONFIG_VERSIONA=`grep ^PACKAGE_VERSION= $CONFIG_FILE | cut -f2 -d= |
cut -f2 -d\'`
CONFIG_VERSION=$CONFIG_MAJOR.$CONFIG_MINOR
test "$CONFIG_PATCH" = 0 || CONFIG_VERSION=$CONFIG_VERSION.$CONFIG_PATCH
test "$CONFIG_VERSION"  = "$VERSION" || echo autoconf version number mismatch
test "$CONFIG_VERSIONA" = "$VERSION" || echo autoconf version string mismatch

cd $TEMP/relx/GeographicLib-$VERSION
(
    echo Files with trailing spaces:
    find . -type f | egrep -v 'config\.guess|Makefile\.in|\.m4|\.png|\.gif|\.pdf' |
	while read f; do
	    tr -d '\r' < $f | grep ' $' > /dev/null && echo $f || true
	done
    echo
    echo Files with tabs:
    find . -type f |
	egrep -v '[Mm]akefile|\.html|\.vcproj|\.sln|\.m4|\.png|\.gif|\.pdf|\.xml' |
	egrep -v '\.sh|depcomp|install-sh|/config\.|configure$|compile|missing' |
	egrep -v 'js/samples/geod-.*\.html' |
	xargs grep -l  '	' || true
    echo
    echo Files with multiple newlines:
    find . -type f |
	egrep -v \
	   '/Makefile\.in|\.1\.html|\.png|\.gif|\.pdf|/ltmain|/config|\.m4|Settings' |
	egrep -v '(Resources|Settings)\.Designer\.cs' |
	while read f; do
	    tr 'X\n' 'xX' < $f | grep XXX > /dev/null && echo $f || true
	done
    echo
    echo Files with no newline at end:
    find . -type f |
	egrep -v '\.png|\.gif|\.pdf' |
	while read f; do
	    n=`tail -1 $f | wc -l`; test $n -eq 0 && echo $f || true
	done
    echo
    echo Files with extra newlines at end:
    find . -type f |
	egrep -v '/configure|/ltmain.sh|\.png|\.gif|\.pdf|\.1\.html' |
	while read f; do
	    n=`tail -1 $f | wc -w`; test $n -eq 0 && echo $f || true
	done
    echo
    echo JS files with bad comment ends:
    find js -type f -name '*.js' | xargs grep -l '\*\*/' || true
    echo
) > $TEMP/badfiles.txt
cat $TEMP/badfiles.txt
cat > $TEMP/tasks.txt <<EOF
# deploy documentation
test -d $WEBDIST/htdocs/$VERSION-pre &&
rm -rf $WEBDIST/htdocs/$VERSION &&
mv $WEBDIST/htdocs/$VERSION{-pre,} &&
make -C $DEVELSOURCE -f makefile-admin distrib-doc

rm $WEBDIST/htdocs/html &&
ln -s $VERSION $WEBDIST/htdocs/html &&
make -C $DEVELSOURCE -f makefile-admin distrib-doc

# deploy release packages
chmod 755 $DEVELSOURCE/GeographicLib-$VERSION-win{32,64}.exe
chmod 644 $DEVELSOURCE/GeographicLib-$VERSION{.tar.gz,.zip}
mv $DEVELSOURCE/GeographicLib-$VERSION{.tar.gz,.zip,-win{32,64}.exe} $DEVELSOURCE/distrib
make -C $DEVELSOURCE -f makefile-admin distrib-files

# install built version
sudo make -C $TEMP/relc/GeographicLib-$VERSION/BUILD-system install

# python release -- authentication via ~/.pypirc
# New method in
# https://packaging.python.org/tutorials/packaging-projects/#uploading-your-project-to-pypi
cd $TEMP/gita/geographiclib/python
python3 setup.py sdist bdist_wheel
python3 -m twine upload dist/*
# installs in /usr/local/lib/python3.7/site-packages/geographiclib
cd; sudo python3 -m pip install --upgrade geographiclib

# java release -- authentication via ~/.m2/settings.xml; this gets signed too
# (multiple ~4 times!).
cd $TEMP/gita/geographiclib/java
mvn clean deploy -P release

# javascript release
# authenticate via .npmrc; _auth value is
#   echo -n cffk:PASSWORD | openssl base64
# decode with
#   echo AUTHSTRING | openssl base64 -d
cd $TEMP/gita/geographiclib/BUILD/js && npm publish geographiclib
make -C $DEVELSOURCE -f makefile-admin distrib-js
make -C $DEVELSOURCE -f makefile-admin install-js
# also update devel branch of node-geographiclib from ??
# git@github.com:yurijmikhalevich/node-geographiclib.git
$TEMP/gita/geographiclib/BUILD/js/geographiclib

# matlab toolbox
chmod 644 $DEVELSOURCE/geographiclib_toolbox_$VERSION.*
mv $DEVELSOURCE/geographiclib_toolbox_$VERSION.* $DEVELSOURCE/matlab-distrib

# commit and tag release branch
cd $TEMP/gitr/geographiclib
git add -A
git commit -m "Version $VERSION ($DATE)"
git tag -m "Version $VERSION ($DATE)" r$VERSION
git push
git push --tags

# tag master branch
cd $DEVELSOURCE
git checkout master
git merge --no-ff $BRANCH -m "Merge from devel for version $VERSION"
git tag -m "Version $VERSION ($DATE)" v$VERSION
git push --all
git push --tags

# Also to do
# post release notices
# set default download files
# make -f makefile-admin distrib-{cgi,html}
# update home brew
#   dir = /usr/local/Homebrew/Library/Taps/homebrew/homebrew-core
#   branch = geographiclib/$VERSION
#   file = Formula/geographiclib.rb
#   brew upgrade geographiclib
#   commit message = geographiclib $VERSION
# upload matlab packages
# update binaries for cgi applications
# trigger build on build-open
EOF
echo cat $TEMP/tasks.txt
cat $TEMP/tasks.txt
