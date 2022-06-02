#!/bin/bash

source ./build-opts.sh

# Edit to suite your system
PREFIX=$GKYLSOFT/OpenBLAS-0.3.15
# Location where dependency sources will be downloaded
DEP_SOURCES=$HOME/gkylsoft/dep_src/

mkdir -p $DEP_SOURCES
cd $DEP_SOURCES

if [ "$DOWNLOAD_PKGS" = "yes" ]
then
    echo "Downloading OpenBLAS .."
    # delete old checkout and builds
    rm -rf OpenBLAS-*
    curl -L https://github.com/xianyi/OpenBLAS/releases/download/v0.3.15/OpenBLAS-0.3.15.tar.gz > OpenBLAS-0.3.15.tar.gz
fi

if [ "$BUILD_PKGS" = "yes" ]
then
    echo "Building OpenBLAS .."
    gunzip -f OpenBLAS-0.3.15.tar.gz
    tar xvf OpenBLAS-0.3.15.tar
    cd OpenBLAS-0.3.15
    make USE_OPENMP=0 NUM_THREADS=1 FC=$FC -j
    make USE_OPENMP=0 NUM_THREADS=1 install PREFIX=$PREFIX -j

    # soft-link 
    ln -sfn $PREFIX $GKYLSOFT/OpenBLAS
fi


