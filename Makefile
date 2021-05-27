
# You can set compiler to use as:
#
# make CC=mpicc 
#

CFLAGS = -O3 -g -Iminus -Izero -Iapps -Ikernels/basis -Ikernels/maxwell -Ikernels/vlasov
PREFIX = ${HOME}/gkylsoft

# Header dependencies
headers = $(wildcard minus/*.h) $(wildcard zero/*.h) $(wildcard apps/*.h) $(wildcard kernels/*/*.h)

# Object files to compile in library
libobjs = $(patsubst %.c,%.o,$(wildcard minus/*.c)) \
	$(patsubst %.c,%.o,$(wildcard zero/*.c)) \
	$(patsubst %.c,%.o,$(wildcard apps/*.c)) \
	$(patsubst %.c,%.o,$(wildcard kernels/*/*.c))

# Make targets: libraties, regression tests and unit tests
all: build/libgkylzero.a \
	$(patsubst %.c,build/%,$(wildcard regression/app_*.c)) build/regression/twostream.ini \
	$(patsubst %.c,build/%,$(wildcard unit/ctest_*.c))

# Library archive
build/libgkylzero.a: ${libobjs} ${headers}
	ar -crs build/libgkylzero.a ${libobjs}

# Regression tests
build/regression/twostream.ini: regression/twostream.ini
	cp regression/twostream.ini build/regression/twostream.ini

build/regression/%: regression/%.c build/libgkylzero.a
	${CC} ${CFLAGS} -o $@ $< -I. -Lbuild -lgkylzero -lm -lpthread 

# Unit tests

build/unit/%: unit/%.c build/libgkylzero.a
	${CC} ${CFLAGS} -o $@ $< -I. -Lbuild -lgkylzero -lm -lpthread 

# Run unit tests
check: $(patsubst %.c,build/%,$(wildcard unit/ctest_*.c))
	./build/unit/ctest_alloc
	./build/unit/ctest_array
	./build/unit/ctest_basis
	./build/unit/ctest_block_topo
	./build/unit/ctest_fv_proj
	./build/unit/ctest_gauss_quad
	./build/unit/ctest_proj_on_basis
	./build/unit/ctest_range
	./build/unit/ctest_rect_apply_bc
	./build/unit/ctest_rect_decomp
	./build/unit/ctest_rect_grid
	./build/unit/ctest_ref_count
	./build/unit/ctest_update_fsm
	./build/unit/ctest_wv_euler
	./build/unit/ctest_wv_iso_euler
	./build/unit/ctest_wv_maxwell
	./build/unit/ctest_wv_ten_moment

install: all
	mkdir -p ${PREFIX}/gkylzero/include
	mkdir -p ${PREFIX}/gkylzero/lib
	mkdir -p ${PREFIX}/gkylzero/bin
	mkdir -p ${PREFIX}/gkylzero/share
	cp ${headers} ${PREFIX}/gkylzero/include
	cp ${headers} build/headers
	rm -f ${PREFIX}/gkylzero/include/gkylzero.h
	@echo "#pragma once\n" > ${PREFIX}/gkylzero/include/gkylzero.h
	@echo "#define TEST_NO_MAIN" >> ${PREFIX}/gkylzero/include/gkylzero.h
	@echo $(patsubst build/headers/%.h,"#include <%.h>\n",$(wildcard build/headers/*.h)) >> ${PREFIX}/gkylzero/include/gkylzero.h
	cp -f build/libgkylzero.a ${PREFIX}/gkylzero/lib
	cp -f build/Makefile.sample ${PREFIX}/gkylzero/share/Makefile
	cp -f regression/app_twostream.c ${PREFIX}/gkylzero/share/app_twostream.c
	cp -f regression/twostream.ini ${PREFIX}/gkylzero/share/twostream.ini
	cp -f build/regression/app_vlasov_kerntm ${PREFIX}/gkylzero/bin/

clean:
	rm -rf build/libgkylzero.a build/regression/twostream.ini build/headers/*.h ${libobjs} build/regression/app_* build/unit/ctest_*

