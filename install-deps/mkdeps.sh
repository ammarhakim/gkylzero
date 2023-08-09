#!/bin/sh

# Defaults
PREFIX=$HOME/g0g2gpu/gkylsoft-amd

# default build options
CC=gcc
CXX=g++
FC=gfortran
MPICC=$PREFIX/openmpi/bin/mpicc
MPICXX=$PREFIX/openmpi/bin/mpicxx

# by default, do not build anything
BUILD_OPENBLAS=no
BUILD_SUPERLU=no
BUILD_SUPERLU_DIST=no
BUILD_SUPERLU_DIST_GPU=no
BUILD_OPENMPI=no

# by default, download as well as build packages
DOWNLOAD_PKGS=yes
BUILD_PKGS=yes

# ----------------------------------------------------------------------------
# Function definitions
# ----------------------------------------------------------------------------

show_help() {
cat <<EOF

./mkdeps.sh CC=cc CXX=cxx FC=gfortran

Build GkylZero dependencies

CC 
CXX                         C and C++ compilers to use
FC                          Fortran compiler to use (only gfortran is supported)
MPICC                       
MPICXX                      MPI C and C++ compilers to use

-h
--help                      This help.

--download                  [yes] Download packages?
--build                     [yes] Build packages?

--prefix=DIR                Prefix where dependencies should be installed.
                            Default is $HOME/gkylsoft

The following flags specify the libraries to build.

--build-openblas            [no] Should we build OpenBLAS?
--build-superlu             [no] Should we build SuperLU (serial)
--build-superlu_dist        [no] Should we build SuperLU (parallel)
--enable-superlu_gpu        [no] Build GPUs lib for SuperLU (needs --build-superlu_dist=yes)
--build-openmpi             [no] Should we build OpenMPI?

EOF
}

# Helper functions

find_program() {
   prog=`command -v "$1" 2>/dev/null`
   if [ -n "$prog" ]
   then
      dirname "$prog"
   fi
}

die() {
   echo "$*"
   echo
   echo "Dependency builds failed."
   echo
   exit 1
}

# ----------------------------------------------------------------------------
# MAIN PROGRAM
# ----------------------------------------------------------------------------

# Parse options

while [ -n "$1" ]
do
   value="`echo $1 | sed 's/[^=]*.\(.*\)/\1/'`"
   key="`echo $1 | sed 's/=.*//'`"
   if `echo "$value" | grep "~" >/dev/null 2>/dev/null`
   then
      echo
      echo '*WARNING*: the "~" sign is not expanded in flags.'
      echo 'If you mean the home directory, use $HOME instead.'
      echo
   fi
   case "$key" in
   -h)
      show_help
      exit 0
      ;;
   --help)
      show_help
      exit 0
      ;;
   --build)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_PKGS="$value"
      ;;      
   --download)
      [ -n "$value" ] || die "Missing value in flag $key."
      DOWNLOAD_PKGS="$value"
      ;;
   CC)
      [ -n "$value" ] || die "Missing value in flag $key."
      CC="$value"
      ;;
   CXX)
      [ -n "$value" ] || die "Missing value in flag $key."
      CXX="$value"
      ;;
   FC)
      [ -n "$value" ] || die "Missing value in flag $key."
      FC="$value"
      ;;   
   --prefix)
      [ -n "$value" ] || die "Missing value in flag $key."
      PREFIX="$value"
      ;;
   --build-openblas)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_OPENBLAS="$value"
      ;;
   --build-superlu)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_SUPERLU="$value"
      ;;
   --build-superlu_dist)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_SUPERLU_DIST="$value"
      ;;
   --enable-superlu_gpu)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_SUPERLU_DIST_GPU="$value"
      ;;
   --build-openmpi)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_OPENMPI="$value"
      ;;   
   *)
      die "Error: Unknown flag: $1"
      ;;
   esac
   shift
done

CMAKE_SUPERLU_DIST_GPU=OFF
# Set package options
if [ "$BUILD_SUPERLU_DIST_GPU" = "yes" ]
then
    CMAKE_SUPERLU_DIST_GPU=ON
fi

# Write out build options for scripts to use
cat <<EOF1 > build-opts.sh
#!/bin/bash

# Generated automatically! Do not edit

# Download/Build options
DOWNLOAD_PKGS=$DOWNLOAD_PKGS
BUILD_PKGS=$BUILD_PKGS

# Installation directory
GKYLSOFT=$PREFIX
# Various compilers
CC=$CC
CXX=$CXX
MPICC=$MPICC
MPICXX=$MPICXX
FC=gfortran

# Package options
CMAKE_SUPERLU_DIST_GPU=$CMAKE_SUPERLU_DIST_GPU

EOF1

build_openblas() {
    if [ "$BUILD_OPENBLAS" = "yes" ]
    then    
	echo "Building OpenBLAS"
	./build-openblas.sh
    fi
}

build_superlu() {
    if [ "$BUILD_SUPERLU" = "yes" ]
    then    
	echo "Building SUPERLU"
	./build-superlu.sh 
    fi
}

build_superlu_dist() {
    if [ "$BUILD_SUPERLU_DIST" = "yes" ]
    then    
	echo "Building SUPERLU Parallel"
	./build-parmetis.sh
	./build-superlu_dist.sh 
    fi
}

build_openmpi() {
    if [ "$BUILD_OPENMPI" = "yes" ]
    then    
	echo "Building OpenMPI"
	./build-openmpi.sh
    fi
}

echo "Installations will be in  $PREFIX"

build_openmpi
build_openblas
build_superlu
build_superlu_dist
