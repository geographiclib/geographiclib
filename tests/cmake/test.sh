#! /bin/sh
for x in net; do
    while read package static cmakeversion stat; do
	builddir=/tmp/test${x}geog-$package-$static-$cmakeversion
	rm -rf $builddir
	mkdir $builddir
	cd $builddir
	cmake=/tmp/cmake-$cmakeversion-win32-x86/bin/cmake
	echo ==================== test${x}geog-$package-$static-$cmakeversion =============
	$cmake -G "Visual Studio 10" -D GeographicLib_USE_STATIC_LIBS=$static -D CMAKE_PREFIX_PATH=c:/tmp/geoginst-$package u:/geographiclib/tests/cmake/test${x}geog
	test -z "$x" &&
	$cmake --build . --config Debug   --target ALL_BUILD
	$cmake --build . --config Release --target ALL_BUILD
	echo XXXXXXXXXXX
	test -z "$x" &&
	./Debug/geodesictest
	./Release/geodesictest
    done <<EOF
2.8.4-STATIC ON  2.8.7
2.8.4-STATIC ON  2.8.9
2.8.4-STATIC ON  2.8.5
2.8.4-STATIC ON  2.8.6
2.8.4-STATIC ON  2.8.8
2.8.4-STATIC ON  2.8.10
2.8.4-STATIC ON  2.8.4
EOF
done
exit

2.8.4-STATIC ON  2.8.11
2.8.4-SHARED OFF 2.8.11
2.8.4-BOTH   ON  2.8.11
2.8.4-BOTH   OFF 2.8.11
2.8.11-STATIC ON  2.8.11
2.8.11-SHARED OFF 2.8.11
2.8.11-BOTH   ON  2.8.11
2.8.11-BOTH   OFF 2.8.11
