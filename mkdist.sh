#!/bin/bash

# This script creates a directory in which all GkylZero code is copied
# and a Makefile to compile them is written. This allows distribution
# of the code as a single compressed directory that can be built with
# a simple make command (hence, not requiring waf or python).

distname=gkylzero-dist

# create directory 
rm -rf $distname
mkdir $distname

# copy included dependency
cp minus/*.h $distname
cp minus/*.c $distname
# copy library code
cp zero/*.h $distname
cp zero/*.c $distname
# copy kernels
cp kernels/*/*.h $distname
cp kernels/*/*.c $distname
# copy app code
cp apps/*.h $distname
cp apps/*.c $distname

cd $distname

# get list of header files
headers=$(ls *.h | tr "\n" " ")
# list of C code for inclusion in library
sources=$(ls *.c | sed 's/\.c/\.o/' | tr "\n" " " )

# back up to root dir to copy regression and unit tests (this is so
# that the regression, unit objects are not archived into the gkylzero
# library)
cd ..
# copy regression tests
cp regression/*.c $distname
# copy unit tests
cp unit/*.c $distname
cd $distname

# list of regression tests
regtests=$(ls app_*.c | sed 's/\.c//' | tr "\n" " " )
# list of unit tests
unittests=$(ls ctest_*.c | sed 's/\.c//' | tr "\n" " " )

# create app targets for insertion in generated Makefile
regtargets=""
for rt in $regtests
do
    regtargets+="$rt: $rt.c libgkylzero.a\n"
    regtargets+="\t\${CC} \${CFLAGS} $rt.c -I. -L. -lgkylzero -lm -lpthread -o $rt\n\n"
done

# create unit-test targets for insertion in generated Makefile
unittargets=""
for ut in $unittests
do
    unittargets+="$ut: $ut.c libgkylzero.a\n"
    unittargets+="\t\${CC} \${CFLAGS} $ut.c -I. -L. -lgkylzero -lm -lpthread -o $ut\n\n"
done

# create Makefile
cat <<EOF > Makefile

# Amalgamated Makefile, generated automatically. Set compiler or
# modify as needed. For example, to set compiler do:
#
# make CC=mpicc 
#

CFLAGS = -O3 -g -I.
PREFIX = \${HOME}/gkylsoft

headers = $headers

libobjs = $sources

all: libgkylzero.a ${regtests} ${unittests}

libgkylzero.a: \${libobjs}
	ar -crs libgkylzero.a \${libobjs}

install: libgkylzero.a
	 mkdir -p \${PREFIX}/gkylzero/include
	 mkdir -p \${PREFIX}/gkylzero/lib
	 mkdir -p \${PREFIX}/gkylzero/bin
	 mkdir -p \${PREFIX}/gkylzero/share
	 cp -f ${headers} \${PREFIX}/gkylzero/include
	 cp -f libgkylzero.a \${PREFIX}/gkylzero/lib
	 cp -f 000version.txt \${PREFIX}/gkylzero
	 cp -f Makefile.sample \${PREFIX}/gkylzero/share/Makefile
	 cp -f app_twostream.c \${PREFIX}/gkylzero/share/app_twostream.c
	 cp -f app_vlasov_kerntm \${PREFIX}/gkylzero/bin/

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

 make PREFIX=mylocation install

The code will be installed in mylocation/gkylzero.

To uninstall manually delete the gkylzero directory wherever you
installed it. There is no explicit uninstall target.

Note that GkylZero is meant to be a library and not a stand-alone
code. To write an application that uses GkylZero a sample app and
Makefile will be installed in the "share" directory of the
installation. Also see app_* files in this directory for example
applications.

Generated on $(date) by $(whoami).

EOF1

# write version information
cat <<EOF2 > 000version.txt
git changeset: $(git describe --abbrev=12 --always --dirty=+)
Archive generated on $(date) by $(whoami).
EOF2

# write sample Makefile for user projects
cat <<EOF3 > Makefile.sample

# Sample Makefile to use installed gkylzero library: copy and modify
# for your needs

CFLAGS = -O3 -g -I.
PREFIX = \${HOME}/gkylsoft

G0_INC_DIR = \${PREFIX}/gkylzero/include
G0_LIB_DIR = \${PREFIX}/gkylzero/lib
G0_LIB_FLAGS = -lm -lgkylzero

all: app_twostream

app_twostream: app_twostream.c
	 \${CC} \${CFLAGS} -I\${G0_INC_DIR} -L\${G0_LIB_DIR} \${G0_LIB_FLAGS} app_twostream.c -o app_twostream

clean:
	rm -rf app_twostream app_twostream.dSYM

EOF3

# up a directory and tar.gz archive
cd ..
tar -zcf ${distname}.tar.gz ${distname}
