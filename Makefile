#
# You can set parameters on the command line:
#
# make CC=mpicc
#
# Or run the configure script to set various parameters. Usually
# defaults are all you need, specially if the dependencies are in
# $HOME/gkylsoft and you are using standard compilers (not building on
# GPUs)

ARCH_FLAGS = -march=native -ffast-math
# Warning flags: -Wall -Wno-unused-variable -Wno-unused-function -Wno-missing-braces
CFLAGS = -O3 -g -ffast-math 
LDFLAGS =
KERN_INCLUDES = -Ikernels/basis -Ikernels/maxwell -Ikernels/vlasov -Ikernels/bin_op -Ikernels/lbo -Ikernels/fem

# Install prefix
PREFIX = ${HOME}/gkylsoft

# Default lapack include and libraries: we prefer linking to static library
LAPACK_INC = ${HOME}/gkylsoft/OpenBLAS/include
LAPACK_LIB = ${HOME}/gkylsoft/OpenBLAS/lib/libopenblas.a

# SuperLU includes and librararies
SUPERLU_INC = ${HOME}/gkylsoft/superlu/include
SUPERLU_LIB = ${HOME}/gkylsoft/superlu/lib/libsuperlu.a

# Include config.mak file (if it exists) to overide defaults above
-include config.mak

# determine OS we are running on
UNAME = $(shell uname)

# On OSX we should use Accelerate framework
ifeq ($(UNAME), Darwin)
	LAPACK_INC = zero # dummy
	LAPACK_LIB = -framework Accelerate
	CFLAGS += -DGKYL_USING_FRAMEWORK_ACCELERATE
endif

INCLUDES = -Iminus -Iminus/STC/include -Izero -Iapps -Iregression ${KERN_INCLUDES} -I${LAPACK_INC} -I${SUPERLU_INC}

NVCC = 
USING_NVCC =
NVCC_FLAGS = 
CUDA_LIBS = 
ifeq ($(CC), nvcc)
       CFLAGS = -O3 -g --forward-unknown-to-host-compiler
       USING_NVCC = yes
       NVCC_FLAGS = -x cu -dc -arch=sm_70 --compiler-options="-fPIC" 
       LDFLAGS += -arch=sm_70
       CUDA_LIBS = -lcublas
endif

%.o : %.cu
	${CC} -c $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -o $@ $<

%.o : %.c
	${CC} -c $(CFLAGS) $(ARCH_FLAGS) $(INCLUDES) -o $@ $< 

# Header dependencies
headers = $(wildcard minus/*.h) $(wildcard zero/*.h) $(wildcard apps/*.h) $(wildcard kernels/*/*.h)

# Object files to compile in library
libobjs = $(patsubst %.c,%.o,$(wildcard minus/*.c)) \
	$(patsubst %.c,%.o,$(wildcard zero/*.c)) \
	$(patsubst %.c,%.o,$(wildcard apps/*.c)) \
	$(patsubst %.c,%.o,$(wildcard kernels/*/*.c))

ifdef USING_NVCC

# Object constructed from cu files to include in library
libobjs += $(patsubst %.cu,%.o,$(wildcard zero/*.cu))

# Unfortunately, due to the limitations of the NVCC compiler to treat
# device code in C files, we need to force compile the kernel code
# using the -x cu flag

kernels/bin_op/%.o : kernels/bin_op/%.c
	${CC} -c $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -o $@ $<

kernels/lbo/%.o : kernels/lbo/%.c
	${CC} -c $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -o $@ $<

kernels/maxwell/%.o : kernels/maxwell/%.c
	${CC} -c $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -o $@ $<

kernels/vlasov/%.o : kernels/vlasov/%.c
	${CC} -c $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -o $@ $<

kernels/basis/%.o : kernels/basis/%.c
	${CC} -c $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -o $@ $<
endif

# Make targets: libraries, regression tests and unit tests
all: build/libgkylzero.a \
	$(patsubst %.c,build/%,$(wildcard regression/rt_*.c)) build/regression/twostream.ini \
	$(patsubst %.c,build/%,$(wildcard unit/ctest_*.c))

# Library archive
build/libgkylzero.a: ${libobjs} ${headers}
	ar -crs build/libgkylzero.a ${libobjs}

# Regression tests
build/regression/twostream.ini: regression/twostream.ini
	cp regression/twostream.ini build/regression/twostream.ini

build/regression/%: regression/%.c build/libgkylzero.a regression/rt_arg_parse.h
	${CC} ${CFLAGS} ${LDFLAGS} -o $@ $< -I. $(INCLUDES) -Lbuild -lgkylzero ${SUPERLU_LIB} ${LAPACK_LIB} ${CUDA_LIBS} -lm -lpthread 

# Unit tests
build/unit/%: unit/%.c build/libgkylzero.a
	${CC} ${CFLAGS} ${LDFLAGS} -o $@ $< -I. $(INCLUDES) -Lbuild -lgkylzero ${SUPERLU_LIB} ${LAPACK_LIB} ${CUDA_LIBS} -lm -lpthread


ifdef USING_NVCC
# unit tests needing CUDA kernels

build/unit/ctest_range: unit/ctest_range.o unit/ctest_range_cu.o build/libgkylzero.a
	${CC} ${LDFLAGS} unit/ctest_range.o unit/ctest_range_cu.o -o build/unit/ctest_range -Lbuild -lgkylzero ${SUPERLU_LIB} ${LAPACK_LIB} ${CUDA_LIBS} -lm -lpthread

build/unit/ctest_rect_grid: unit/ctest_rect_grid.o unit/ctest_rect_grid_cu.o build/libgkylzero.a
	${CC} ${LDFLAGS} unit/ctest_rect_grid.o unit/ctest_rect_grid_cu.o -o build/unit/ctest_rect_grid -Lbuild -lgkylzero ${SUPERLU_LIB} ${LAPACK_LIB} ${CUDA_LIBS} -lm -lpthread

build/unit/ctest_array: unit/ctest_array.o unit/ctest_array_cu.o build/libgkylzero.a
	${CC} ${LDFLAGS} unit/ctest_array.o unit/ctest_array_cu.o -o build/unit/ctest_array -Lbuild -lgkylzero ${SUPERLU_LIB} ${LAPACK_LIB} ${CUDA_LIBS} -lm -lpthread

build/unit/ctest_dg_bin_op: unit/ctest_dg_bin_op.o unit/ctest_dg_bin_op_cu.o build/libgkylzero.a
	${CC} ${LDFLAGS} unit/ctest_dg_bin_op.o unit/ctest_dg_bin_op_cu.o -o build/unit/ctest_dg_bin_op -Lbuild -lgkylzero ${SUPERLU_LIB} ${LAPACK_LIB} ${CUDA_LIBS} -lm -lpthread

build/unit/ctest_dg_maxwell: unit/ctest_dg_maxwell.o unit/ctest_dg_maxwell_cu.o build/libgkylzero.a
	${CC} ${LDFLAGS} unit/ctest_dg_maxwell.o unit/ctest_dg_maxwell_cu.o -o build/unit/ctest_dg_maxwell -Lbuild -lgkylzero ${SUPERLU_LIB} ${LAPACK_LIB} ${CUDA_LIBS} -lm -lpthread

build/unit/ctest_dg_vlasov: unit/ctest_dg_vlasov.o unit/ctest_dg_vlasov_cu.o build/libgkylzero.a
	${CC} ${LDFLAGS} unit/ctest_dg_vlasov.o unit/ctest_dg_vlasov_cu.o -o build/unit/ctest_dg_vlasov -Lbuild -lgkylzero ${SUPERLU_LIB} ${LAPACK_LIB} ${CUDA_LIBS} -lm -lpthread

build/unit/ctest_mom_vlasov: unit/ctest_mom_vlasov.o unit/ctest_mom_vlasov_cu.o build/libgkylzero.a
	${CC} ${LDFLAGS} unit/ctest_mom_vlasov.o unit/ctest_mom_vlasov_cu.o -o build/unit/ctest_mom_vlasov -Lbuild -lgkylzero ${SUPERLU_LIB} ${LAPACK_LIB} ${CUDA_LIBS} -lm -lpthread

build/unit/ctest_hyper_dg: unit/ctest_hyper_dg.o unit/ctest_hyper_dg_cu.o build/libgkylzero.a
	${CC} ${LDFLAGS} unit/ctest_hyper_dg.o unit/ctest_hyper_dg_cu.o -o build/unit/ctest_hyper_dg -Lbuild -lgkylzero ${SUPERLU_LIB} ${LAPACK_LIB} ${CUDA_LIBS} -lm -lpthread

endif

.PHONY: check clean partclean install

# Run unit tests
check: $(patsubst %.c,build/%,$(wildcard unit/ctest_*.c))
	./build/unit/ctest_alloc
	./build/unit/ctest_array
	./build/unit/ctest_array_reduce
	./build/unit/ctest_basis
	./build/unit/ctest_bcorr
	./build/unit/ctest_block_topo
	./build/unit/ctest_dg_bin_ops
	./build/unit/ctest_dg_lbo_vlasov
	./build/unit/ctest_dg_maxwell
	./build/unit/ctest_dg_vlasov
	./build/unit/ctest_dynvec
	./build/unit/ctest_fv_proj
	./build/unit/ctest_gauss_quad
	./build/unit/ctest_hyper_dg
	./build/unit/ctest_mat
	./build/unit/ctest_mat_triples
	./build/unit/ctest_mom_calc
	./build/unit/ctest_mom_vlasov
	./build/unit/ctest_proj_maxwellian_on_basis
	./build/unit/ctest_proj_on_basis
	./build/unit/ctest_range
	./build/unit/ctest_rect_apply_bc
	./build/unit/ctest_rect_decomp
	./build/unit/ctest_rect_grid
	./build/unit/ctest_ref_count
	./build/unit/ctest_superlu
	./build/unit/ctest_update_fsm
	./build/unit/ctest_wave_geom
	./build/unit/ctest_wv_euler
	./build/unit/ctest_wv_iso_euler
	./build/unit/ctest_wv_maxwell
	./build/unit/ctest_wv_mhd
	./build/unit/ctest_wv_sr_euler
	./build/unit/ctest_wv_ten_moment

install: all
	mkdir -p ${PREFIX}/gkylzero/include
	mkdir -p ${PREFIX}/gkylzero/lib
	mkdir -p ${PREFIX}/gkylzero/bin
	mkdir -p ${PREFIX}/gkylzero/share
	cp ${headers} ${PREFIX}/gkylzero/include
	cp -f build/libgkylzero.a ${PREFIX}/gkylzero/lib
	cp -f build/Makefile.sample ${PREFIX}/gkylzero/share/Makefile
	cp -f regression/rt_arg_parse.h ${PREFIX}/gkylzero/share/rt_arg_parse.h
	cp -f regression/rt_twostream.c ${PREFIX}/gkylzero/share/rt_twostream.c
	cp -f regression/twostream.ini ${PREFIX}/gkylzero/share/twostream.ini
	cp -f build/regression/rt_vlasov_kerntm ${PREFIX}/gkylzero/bin/

clean:

	rm -rf build/libgkylzero.a build/regression/twostream.ini */*.o kernels/*/*.o build/regression/rt_* build/unit/ctest_*

# partclean does not delete the kernel object files as they do not
# always need to be rebuilt and are very expensive to build when using
# nvcc
partclean:
	rm -rf build/libgkylzero.a app/*.o zero/*.o unit/*.o build/regression/rt_* build/unit/ctest_*
