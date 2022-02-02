#!/bin/sh

# Defaults
PREFIX=$HOME/gkylsoft

# default build options
CC=gcc
CXX=g++
FC=gfortran

# by default, do not build anything
BUILD_OPENBLAS=no
BUILD_SUPERLU=no

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

-h
--help                      This help.

--download                  [yes] Download packages?
--build                     [yes] Build packages?

--prefix=DIR                Prefix where dependencies should be installed.
                            Default is $HOME/gkylsoft

The following flags specify the libraries to build.

--build-openblas            [no] Should we build OpenBLAS?
--build-superlu             [no] Should we build SuperLU (serial)

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
   *)
      die "Error: Unknown flag: $1"
      ;;
   esac
   shift
done

# Write out build options for scripts to use
cat <<EOF1 > build-opts.sh
# Generated automatically! Do not edit

# Download/Build options
DOWNLOAD_PKGS=$DOWNLOAD_PKGS
BUILD_PKGS=$BUILD_PKGS

# Installation directory
GKYLSOFT=$PREFIX
# Various compilers
CC=$CC
CXX=$CXX
FC=gfortran

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

echo "Installations will be in  $PREFIX"

build_openblas
build_superlu

