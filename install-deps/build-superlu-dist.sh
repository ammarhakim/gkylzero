#!/bin/bash

source ./build-opts.sh

# Edit to suite your system
PREFIX=$GKYLSOFT/superlu-
# Location where dependency sources will be downloaded
DEP_SOURCES=$HOME/gkylsoft/dep_src/superlu_dist-8.1.0

mkdir -p $DEP_SOURCES
cd $DEP_SOURCES

if [ "$DOWNLOAD_PKGS" = "yes" ]
then
    echo "Downloading SuperLU Distributed .."
    # delete old checkout and builds
    rm -rf superlu_dist-*
    curl -L https://github.com/xiaoyeli/superlu_dist/archive/refs/tags/v8.1.0.tar.gz > superlu_dist-8.1.0.tar.gz
fi

if [ "$BUILD_PKGS" = "yes" ]
then
    echo "Building SuperLU .."
    gunzip -f superlu_dist-8.1.0.tar.gz
    tar xvf superlu_dist-8.1.0.tar

    cd superlu_dist-8.1.0
    mkdir build
    cd build

    # configure build
#    cmake .. -DCMAKE_C_FLAGS="-g -O3 -fPIC" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_INSTALL_LIBDIR=lib -Denable_tests=NO -Denable_internal_blaslib=NO -DXSDK_ENABLE_Fortran=NO

    # build and install
#    make -j 32 VERBOSE=1
#    make install

    # soft-link 
#    ln -sfn $PREFIX $GKYLSOFT/superlu
fi
