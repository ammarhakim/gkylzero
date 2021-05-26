
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
all: makeout/libgkylzero.a \
	$(patsubst %.c,makeout/%,$(wildcard regression/app_*.c)) \
	$(patsubst %.c,makeout/%,$(wildcard unit/ctest_*.c))

# Library archive
makeout/libgkylzero.a: ${libobjs}
	ar -crs makeout/libgkylzero.a ${libobjs}

# Regression tests
makeout/regression/twostream.ini: regression/twostream.ini
	cp regression/twostream.ini makeout/regression/twostream.ini

makeout/regression/%: regression/%.c makeout/libgkylzero.a
	${CC} ${CFLAGS} -I. -Lmakeout -lgkylzero -lm -lpthread -o $@ $<

# Unit tests

makeout/unit/%: unit/%.c makeout/libgkylzero.a
	${CC} ${CFLAGS} -I. -Lmakeout -lgkylzero -lm -lpthread -o $@ $<

# Run unit tests
check: $(patsubst %.c,makeout/%,$(wildcard unit/ctest_*.c))
	./makeout/unit/ctest_alloc
	./makeout/unit/ctest_array
	./makeout/unit/ctest_basis
	./makeout/unit/ctest_block_topo
	./makeout/unit/ctest_fv_proj
	./makeout/unit/ctest_gauss_quad
	./makeout/unit/ctest_proj_on_basis
	./makeout/unit/ctest_range
	./makeout/unit/ctest_rect_apply_bc
	./makeout/unit/ctest_rect_decomp
	./makeout/unit/ctest_rect_grid
	./makeout/unit/ctest_ref_count
	./makeout/unit/ctest_update_fsm
	./makeout/unit/ctest_wv_euler
	./makeout/unit/ctest_wv_iso_euler
	./makeout/unit/ctest_wv_maxwell
	./makeout/unit/ctest_wv_ten_moment

clean:
	rm -rf makeout/libgkylzero.a ${libobjs} makeout/regression/app_* makeout/unit/ctest_*

