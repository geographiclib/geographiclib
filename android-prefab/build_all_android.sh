#!/bin/bash

set -e

if [ -z "$ANDROID_NDK" ]; then
  echo "Please set ANDROID_NDK to the Android NDK folder"
  exit 1
fi

# Directories, paths and filenames
BUILD_DIR=build

CMAKE_ARGS="-H..\ \
    -DCMAKE_BUILD_TYPE=Release \
    -DANDROID_STL=c++_shared \
    -DCMAKE_TOOLCHAIN_FILE=${ANDROID_NDK}/build/cmake/android.toolchain.cmake \
    -DANDROID_NDK=${ANDROID_NDK} \
    -DBUILD_SHARED_LIBS=OFF"

rm -rf $BUILD_DIR

function build_geographiclib {

  ABI=$1
  ABI_BUILD_DIR=build/${ABI}

  echo "Building GeographicLib for ${ABI}"

  mkdir -p ${ABI_BUILD_DIR} ${ABI_BUILD_DIR}/${STAGING_DIR}

  cmake -B${ABI_BUILD_DIR} \
        -DANDROID_ABI=${ABI} \
        -DCMAKE_INSTALL_PREFIX=build/staging/${ABI} \
        ${CMAKE_ARGS}

  cmake --build ${ABI_BUILD_DIR} \
        --target install \
        -j
}

build_geographiclib armeabi-v7a
build_geographiclib arm64-v8a
build_geographiclib x86
build_geographiclib x86_64
