#!/bin/bash

source ./build-opts.sh

# Install prefix
PREFIX=$GKYLSOFT/openmpi-4.0.5

# delete old checkout and builds
rm -rf openmpi-4.0.5.tar* openmpi-4.0.5

curl -L https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.5.tar.gz > openmpi-4.0.5.tar.gz
gunzip -c openmpi-4.0.5.tar.gz | tar xf -
cd openmpi-4.0.5
./configure --prefix=$PREFIX --enable-mpi-fortran=none CC=$CC CXX=$CXX
make all install

# soft-link 
ln -sfn $PREFIX $GKYLSOFT/openmpi
