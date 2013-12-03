#! /bin/sh -e

while read v; do
    for t in STATIC SHARED BOTH; do
	(
	    echo ================================ $v-$t ============
	    b=geogbuild-$v-$t
	    inst=geoginst-$v-$t/GeographicLib-1.34
	    rm -rf /tmp/$b /tmp/$inst
	    mkdir /tmp/$b
	    cd /tmp/$b
	    cmake=/tmp/cmake-$v-win32-x86/bin/cmake
	    $cmake -G "Visual Studio 10" -D GEOGRAPHICLIB_LIB_TYPE=$t -D CMAKE_INSTALL_PREFIX=c:/tmp/$inst -D BUILD_NETGEOGRAPHICLIB=ON u:/geographiclib
	    $cmake --build . --config Debug --target ALL_BUILD
	    $cmake --build . --config Debug --target RUN_TESTS
	    $cmake --build . --config Debug --target INSTALL
	    $cmake --build . --config Release --target ALL_BUILD
	    $cmake --build . --config Release --target RUN_TESTS
	    $cmake --build . --config Release --target INSTALL
	)
    done
done <<EOF
2.8.4
2.8.11
EOF

exit

2.8.0
2.8.1
2.8.2
2.8.3
2.8.4
2.8.5
2.8.6
2.8.7
2.8.8
2.8.9
2.8.10
2.8.10.1
2.8.10.2
2.8.11
2.8.11.1
2.8.11.2
2.8.12
2.8.12.1
