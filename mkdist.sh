#!/bin/bash

# This script creates a directory in which all GkylZero code is copied
# and a Makefile to compile them is written. This allows distribution
# of the code as a single compressed directory that can be built with
# a simple make command (hence, not requiring waf or python).

distname=gzero

# create directory 
rm -rf $distname
mkdir $distname

# copy library code
cp zero/*.h $distname
cp zero/*.c $distname
# copy app code
cp miniapp/*.h $distname
cp miniapp/*.c $distname
# copy kernels
cp kernels/*/*.h $distname
cp kernels/*/*.c $distname
# copy regression tests
cp regression/*.c $distname
# copy unit tests
cp unit/acutest.h $distname
cp unit/*.c $distname

cd $distname

# get list of header files
headers=$(ls *.h | tr "\n" " ")
# list of C code
sources=$(ls *.c | sed 's/\.c/\.o/' | tr "\n" " " )
# list of regression tests
regtests=$(ls app_*.c | sed 's/\.c//' | tr "\n" " " )
# list of unit tests
unittests=$(ls ctest_*.c | sed 's/\.c//' | tr "\n" " " )

# create app targets for insertion in generated Makefile
regtargets=""
for rt in $regtests
do
    regtargets+="$rt: $rt.c libgkylzero.a\n"
    regtargets+="\t\${CC} \${CFLAGS} $rt.c -I. -L. -lgkylzero -lm -o $rt\n\n"
done

# create unit-test targets for insertion in generated Makefile
unittargets=""
for ut in $unittests
do
    unittargets+="$ut: $ut.c libgkylzero.a\n"
    unittargets+="\t\${CC} \${CFLAGS} $ut.c -I. -L. -lgkylzero -lm -o $ut\n\n"
done

# create Makefile
cat <<EOF > Makefile

# Amalgamated makefile. Set compiler if needed:
#
# make CC=mpicc 
#

CFLAGS = -O3 -g -I.

headers = $headers

libobjs = $sources

all: libgkylzero.a ${regtests} ${unittests}

libgkylzero.a: \${libobjs}
	ar -crs libgkylzero.a \${libobjs}

$(printf "%b" "$regtargets")

$(printf "%b" "$unittargets")

clean:
	rm -rf libgkylzero.a \${libobjs} ${regtests} ${unittests} *.dSYM *.so 

EOF

# write a README file
cat <<EOF1 > README.md

This directory contains the complete source for GkylZero copied to a
single directory. It also contains an auto-generated
Makefile. Documentation for Gkeyll is available at
https://gkeyll.readthedocs.io.

# Install

You can build the code by simply typing 'make' in this directory. As
GkylZero has minimal external dependencies, this should "just work" on
most systems. However, sometimes you might need to pass options to the
make command. For example, to use Intel compiler do:

 make CC=icc

You can of course edit the Makefile yourself to suite your needs. To
install the code do:

 make install

The default installation location is \$HOME/gkylsoft/gkylzero. You can
change it by specifying the PREFIX variable:

 make PREFIX=/mylocation install

Generated on $(date) by $(whoami).

EOF1

# up a directory and tar.gz archive
cd ..
tar -zcf ${distname}.tar.gz ${distname}
