#!/bin/bash

set -e

# Test that the prefab binary exists
if hash prefab 2>/dev/null; then
  echo "Prefab is installed"
else
  echo "Prefab binary not found. See https://github.com/google/prefab#building-the-command-line-executable for install instructions"
  exit 1
fi

# Get the version string from the source
major=$(sed -n 's/set (PROJECT_VERSION_MAJOR \([0-9]*\))/\1/p' ../CMakeLists.txt)
minor=$(sed -n 's/set (PROJECT_VERSION_MINOR \([0-9]*\))/\1/p' ../CMakeLists.txt)
patch=$(sed -n 's/set (PROJECT_VERSION_PATCH \([0-9]*\))/\1/p' ../CMakeLists.txt)
version="${major}.${minor}.${patch}"

echo "Building libraries for geographiclib version "$version
. build_all_android.sh

rm -rf build/prefab
mkdir -p build/prefab
cp -R prefab/* build/prefab

ABIS=("x86" "x86_64" "arm64-v8a" "armeabi-v7a")

pushd build/prefab

  # Write the version number into the various metadata files
  mv geographiclib-VERSION geographiclib-$version
  mv geographiclib-VERSION.pom geographiclib-$version.pom
  sed -i -e "s/VERSION/${version}/g" geographiclib-$version.pom geographiclib-$version/prefab/prefab.json

  for abi in ${ABIS[@]}
  do
    echo "Copying the ${abi} library"
    # Copy the headers
    cp -R "../staging/${abi}/include" "geographiclib-$version/prefab/modules/geographiclib/libs/android.${abi}/"
    # Copy the libraries
    cp -v "../staging/${abi}/lib/libGeographicLib.a" "geographiclib-${version}/prefab/modules/geographiclib/libs/android.${abi}/"
  done

  # Verify the prefab packages
  for abi in ${ABIS[@]}
  do

    prefab --build-system cmake --platform android --os-version 31 \
        --stl c++_shared --ndk-version 24 --abi ${abi} \
        --output prefab-output-tmp $(pwd)/geographiclib-${version}/prefab

    result=$?;
    if [[ $result == 0 ]]; then
      echo "${abi} package verified"
    else
      echo "${abi} package verification failed"
      exit 1
    fi
  done

  # Zip into an AAR and move into parent dir
  pushd geographiclib-${version}
    zip -r geographiclib-${version}.aar . 2>/dev/null;
    zip -Tv geographiclib-${version}.aar 2>/dev/null;

    # Verify that the aar contents are correct (see output below to verify)
    result=$?;
    if [[ $result == 0 ]]; then
      echo "AAR verified"
    else
      echo "AAR verification failed"
      exit 1
    fi

    mv geographiclib-${version}.aar ..
  popd

  # Zip the .aar and .pom files into a maven package
  zip geographiclib-${version}.zip geographiclib-${version}.* 2>/dev/null;
  zip -Tv geographiclib-${version}.zip 2>/dev/null;

  # Verify that the zip contents are correct (see output below to verify)
  result=$?;
  if [[ $result == 0 ]]; then
    echo "Zip verified"
  else
    echo "Zip verification failed"
    exit 1
  fi

  echo "Prefab zip ready for deployment: ./build/prefab/geographiclib-${version}.zip"
popd
