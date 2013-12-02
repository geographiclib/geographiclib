#! /bin/sh
for x in "" net; do
while read package libtype cmakeversion; do
echo ==================== ${x}geogtest-$package-$libtype-$cmakeversion =============
    builddir=/tmp/${x}geogtest-$package-$libtype-$cmakeversion
    rm -rf $builddir
    mkdir $builddir
    cd $builddir
    case $cmakeversion in
	old ) cmake=/u/temp/cmake-2.8.0/bin/cmake;;
	new ) cmake=cmake;;
    esac
    case $libtype in
	shared ) static=OFF;;
	static ) static=ON;
    esac
    $cmake -G "Visual Studio 10" -D GeographicLib_USE_STATIC_LIBS=$static -D CMAKE_PREFIX_PATH=c:/tmp/$package u:/geographiclib/tests/cmake/test${x}geog
    test -z "$x" && $cmake --build . --config Debug --target ALL_BUILD
    $cmake --build . --config Release --target ALL_BUILD
    test -z "$x" && ./Debug/geodesictest
    ./Release/geodesictest
done <<EOF
geoginst0 shared new
geoginst0 shared old
geoginst0 static new
geoginst0 static old
geoginst1 shared new
geoginst1 shared old
geoginst1 static new
geoginst1 static old
EOF
done
exit
