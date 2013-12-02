#! /bin/sh -e

for v in 0 1; do
    (
	b=geogbuild$v
	inst=geoginst$v/GeographicLib-1.34
	rm -rf /tmp/$b /tmp/$inst
	mkdir /tmp/$b
	cd /tmp/$b
	case $v in
	    0 ) cmake=/u/temp/cmake-2.8.0/bin/cmake;;
	    1 ) cmake=cmake;;
	esac
	$cmake -G "Visual Studio 10" -D GEOGRAPHICLIB_LIB_TYPE=BOTH -D CMAKE_INSTALL_PREFIX=c:/tmp/$inst -D BUILD_NETGEOGRAPHICLIB=ON u:/geographiclib
	$cmake --build . --config Debug --target ALL_BUILD
	$cmake --build . --config Debug --target RUN_TESTS
	$cmake --build . --config Debug --target INSTALL
	$cmake --build . --config Release --target ALL_BUILD
	$cmake --build . --config Release --target RUN_TESTS
	$cmake --build . --config Release --target INSTALL
    )
done
