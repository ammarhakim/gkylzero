#!/bin/bash

source ./build-opts.sh

# Edit to suite your system
PREFIX=$GKYLSOFT/superlu-5.2.2
# Location where dependency sources will be downloaded
DEP_SOURCES=$HOME/gkylsoft/dep_src/

mkdir -p $DEP_SOURCES
cd $DEP_SOURCES

# delete old checkout and builds
rm -rf superlu-*

curl -L https://github.com/xiaoyeli/superlu/archive/refs/tags/v5.2.2.tar.gz > superlu-5.2.2.tar.gz
gunzip -f superlu-5.2.2.tar.gz
tar xvf superlu-5.2.2.tar

cd superlu-5.2.2
mkdir build
cd build

# configure build
cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX -Denable_tests=NO -Denable_internal_blaslib=NO 

# build and install
make -j
make install

# soft-link 
ln -sfn $PREFIX $GKYLSOFT/superlu
