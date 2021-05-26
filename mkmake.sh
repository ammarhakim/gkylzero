#!/bin/bash

# This script creates a Makefile for use in building the code.

# list of library headers and sources
lib_headers=$(ls minus/*.h zero/*.h apps/*.h | tr "\n" " ")
lib_sources=$(ls minus/*.c zero/*.c apps/*.c | sed 's/\.c/\.o/' | tr "\n" " " )

# list of kernel headers and sources
ker_headers=$(ls kernels/*/*.h | tr "\n" " ")
ker_sources=$(ls kernels/*/*.c | sed 's/\.c/\.o/' | tr "\n" " " )

# list of regression and unit tests
regtests=$(ls regression/app_*.c | sed 's/\.c//' | tr "\n" " " )
unittests=$(ls unit/ctest_*.c | sed 's/\.c//' | tr "\n" " " )

# create app targets for insertion in generated Makefile
regtargets=""
for rt in $regtests
do
    regtargets+="makeout/$rt: $rt.c makeout/libgkylzero.a\n"
    regtargets+="\t\${CC} \${CFLAGS} $rt.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/$rt\n\n"
done

# create unit-test targets for insertion in generated Makefile
unittargets=""
for ut in $unittests
do
    unittargets+="makeout/$ut: $ut.c makeout/libgkylzero.a\n"
    unittargets+="\t\${CC} \${CFLAGS} $ut.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/$ut\n\n"
done

# targets list to insert into Makefile
reg_targlist=""
for rt in $regtests
do
    reg_targlist+="makeout/$rt "
done

unit_targlist=""
for ut in $unittests
do
    unit_targlist+="makeout/$ut "
done

# add unit-tests for use in "make check"
checkcmds=""
for ut in $unittests
do
    checkcmds+="\t./makeout/$ut\n"
done

# create Makefile
cat <<EOF > Makefile

# Makefile, generated automatically. Avoid modifying by hand. To
# regenerate this type "./mkmake.sh" in the top-level directory. You
# typically only need to do this if you are adding new code to
# GkeyllZero. You can set compiler to use as:
#
# make CC=mpicc 
#

CFLAGS = -O3 -g -Iminus -Izero -Iapps -Ikernels/basis -Ikernels/maxwell -Ikernels/vlasov
PREFIX = \${HOME}/gkylsoft

headers = ${lib_headers} ${ker_headers}

libobjs = ${lib_sources} ${ker_sources}

all: makeout/libgkylzero.a ${reg_targlist} ${unit_targlist} makeout/regression/twostream.ini

makeout/libgkylzero.a: \${libobjs}
	ar -crs makeout/libgkylzero.a \${libobjs}

makeout/regression/twostream.ini: regression/twostream.ini
	cp regression/twostream.ini makeout/regression/twostream.ini

$(printf "%b" "$regtargets")

$(printf "%b" "$unittargets")

check: ${unit_targlist}
$(printf "%b" "$checkcmds")

clean:
	rm -rf makeout/libgkylzero.a \${libobjs} makeout/regression/app_* makeout/regression/cunit_*

EOF
