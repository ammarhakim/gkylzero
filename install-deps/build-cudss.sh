#!/bin/bash

source ./build-opts.sh

# Edit to suite your system
PREFIX=$GKYLSOFT/cuDSS-0.2.1.3
# Location where dependency sources will be downloaded
DEP_SOURCES=$GKYLSOFT/dep_src/

mkdir -p $DEP_SOURCES
cd $DEP_SOURCES

if [ "$DOWNLOAD_PKGS" = "yes" ]
then
    echo "Downloading cuDSS .."
    # delete old checkout and builds
    rm -rf cuDSS-*
    wget https://developer.download.nvidia.com/compute/cudss/redist/libcudss/linux-x86_64/libcudss-linux-x86_64-0.2.1.3_cuda12-archive.tar.xz
fi

if [ "$BUILD_PKGS" = "yes" ]
then
    echo "Building cuDSS .."
    tar xf libcudss-linux-x86_64-0.2.1.3_cuda12-archive.tar.xz
    mv libcudss-linux-x86_64-0.2.1.3_cuda12-archive $PREFIX

    # soft-link 
    ln -sfn $PREFIX $GKYLSOFT/cuDSS
fi
