#!/bin/bash

source ./build-opts.sh

# Edit to suite your system
PREFIX=$GKYLSOFT/superlu_dist-5.2.2
# Location where dependency sources will be downloaded
DEP_SOURCES=$GKYLSOFT/dep_src/

mkdir -p $DEP_SOURCES
cd $DEP_SOURCES

if [ "$DOWNLOAD_PKGS" = "yes" ]
then
    echo "Downloading SuperLU Distributed .."
    # delete old checkout and builds
    rm -rf superlu_dist-*
    curl -L https://github.com/xiaoyeli/superlu_dist/archive/refs/tags/v8.1.0.tar.gz > superlu_dist-8.1.0.tar.gz
fi

PARMETISLIB=libparmetis.so
OSTYPE=`uname`
if [ "$OSTYPE" = "Darwin" ]
then
    PARMETISLIB=libparmetis.dylib
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
    CC=$MPICC CXX=$MPICXX cmake .. -DCMAKE_C_FLAGS="-g -O3 -fPIC" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_INSTALL_LIBDIR=lib -Denable_tests=NO -Denable_openmp=NO -DXSDK_ENABLE_Fortran=NO -DTPL_ENABLE_PARMETISLIB=YES -DTPL_PARMETIS_LIBRARIES=$GKYLSOFT/parmetis/lib/$PARMETISLIB -DTPL_PARMETIS_INCLUDE_DIRS=$GKYLSOFT/parmetis/include;$GKYLSOFT/metis/include -DTPL_ENABLE_CUDALIB=$CMAKE_SUPERLU_DIST_GPU

    # build and install
    make -j 32 VERBOSE=1
    make install

    # soft-link 
    ln -sfn $PREFIX $GKYLSOFT/superlu_dist
fi
