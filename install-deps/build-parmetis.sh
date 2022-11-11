#!/bin/bash

source ./build-opts.sh

# Edit to suite your system
PREFIX=$GKYLSOFT/parmetis-4.0.3
# Location where dependency sources will be downloaded
DEP_SOURCES=$GKYLSOFT/dep_src/

mkdir -p $DEP_SOURCES
cd $DEP_SOURCES

if [ "$DOWNLOAD_PKGS" = "yes" ]
then
    echo "Downloading ParMetis Distributed .."
    # delete old checkout and builds
    rm -rf superlu_dist-*
    curl -L http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz > parmetis-4.0.3.tar.gz
fi

if [ "$BUILD_PKGS" = "yes" ]
then
    echo "Building ParMetis .."
    gunzip -f parmetis-4.0.3.tar.gz
    tar xvf parmetis-4.0.3.tar

    cd parmetis-4.0.3

    # Build Metis firs
    cd metis
    make config shared=1 prefix=$GKYLSOFT/parmetis-4.0.3
    make -j 32 
    make install

    cd ..
    # configure build
    make config shared=1 prefix=$GKYLSOFT/parmetis-4.0.3 cc=$MPICC cxx=$MPICXX

    # build and install
    make -j 32 
    make install

    # soft-link 
    ln -sfn $PREFIX $GKYLSOFT/parmetis
fi
