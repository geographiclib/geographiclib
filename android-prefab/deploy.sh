#!/bin/bash

if hash mvn 2>/dev/null; then
  echo "maven is installed"
else
  echo "maven binary not found."
  exit 1
fi

# to avoid "gpg: signing failed: Inappropriate ioctl for device" error when run through ssh ...
export GPG_TTY=$(tty)

# Get the version string from the source
major=$(sed -n 's/set (PROJECT_VERSION_MAJOR \([0-9]*\))/\1/p' ../CMakeLists.txt)
minor=$(sed -n 's/set (PROJECT_VERSION_MINOR \([0-9]*\))/\1/p' ../CMakeLists.txt)
patch=$(sed -n 's/set (PROJECT_VERSION_PATCH \([0-9]*\))/\1/p' ../CMakeLists.txt)
version="${major}.${minor}.${patch}"

pushd build/prefab

    echo "Deploy to maven Local repository"
    mvn "install:install-file" \
        "-DpomFile=geographiclib-${version}.pom" \
        "-Dfile=geographiclib-${version}.aar"

#https://central.sonatype.org/pages/manual-staging-bundle-creation-and-deployment.html#manually-deploying-to-ossrh-introduction

    echo "Deploy to sonatype OSSRH maven remote repository"
    mvn "gpg:sign-and-deploy-file" \
        "-Durl=https://oss.sonatype.org/service/local/staging/deploy/maven2/" \
        "-DpomFile=geographiclib-${version}.pom" \
        "-Dfile=geographiclib-${version}.aar" \
        "-DrepositoryId=ossrh" \
        "--settings" "${HOME}/.m2/settings.xml"

popd
