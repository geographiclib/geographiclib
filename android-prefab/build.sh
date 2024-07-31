#!/bin/bash

set -e

#add prefab binary to PATH
PATH=$PATH:$HOME/devel/prefab/cli/build/install/prefab/bin/

#Path to Android NDK
export ANDROID_NDK=$HOME/Android/Sdk/ndk/25.1.8937393/

pushd  "$( dirname -- "$0"; )"

    # create prefab aar package
    . prefab_build.sh

    # deploy package to Local mvm repostory and to OSSRH
    . deploy.sh

popd
