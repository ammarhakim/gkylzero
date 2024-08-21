# -*- makefile-gmake -*-

# Type make help to see help for this Makefile

ARCH_FLAGS ?= -march=native
CUDA_ARCH ?= 70
# Warning flags: -Wall -Wno-unused-variable -Wno-unused-function -Wno-missing-braces
CFLAGS ?= -O3 -g -ffast-math -fPIC -MMD -MP 
LDFLAGS = 
PREFIX ?= ${HOME}/gkylsoft
DEP_PREFIX ?= ${HOME}/gkyl_gpu_dependencies
INSTALL_PREFIX ?= ${PREFIX}

# determine OS we are running on
UNAME = $(shell uname)

# Default lapack include and libraries: we prefer linking to static library
LAPACK_INC = $(DEP_PREFIX)/OpenBLAS/include
LAPACK_LIB_DIR = $(DEP_PREFIX)/OpenBLAS/lib
LAPACK_LIB = -lopenblas

# SuperLU includes and librararies
SUPERLU_INC = $(DEP_PREFIX)/superlu/include
ifeq ($(UNAME_S),Linux)
	SUPERLU_LIB_DIR = $(DEP_PREFIX)/superlu/lib64
	SUPERLU_LIB = $(DEP_PREFIX)/superlu/lib64/libsuperlu.a
else
	SUPERLU_LIB_DIR = $(DEP_PREFIX)/superlu/lib
	SUPERLU_LIB = $(DEP_PREFIX)/superlu/lib/libsuperlu.a
endif

# Include config.mak file (if it exists) to overide defaults above
-include config.mak

# By default, build the "all" target. This builds the G0 shared
# library only. Unit and regression tests are built with explicit
# targets. See "make help"
.DEFAULT_GOAL := all

# CUDA flags
USING_NVCC =
NVCC_FLAGS = 
CUDA_LIBS =
ifeq ($(CC), nvcc)
       USING_NVCC = yes
       CFLAGS = -O3 -g --forward-unknown-to-host-compiler --use_fast_math -ffast-math -MMD -MP -fPIC
       NVCC_FLAGS = -x cu -dc -arch=sm_${CUDA_ARCH} --compiler-options="-fPIC"
       LDFLAGS += -arch=sm_${CUDA_ARCH}
       ifdef CUDAMATH_LIBDIR
              CUDA_LIBS = -L${CUDAMATH_LIBDIR}
       else
              CUDA_LIBS =
       endif
       CUDA_LIBS += -lcublas -lcusparse -lcusolver
endif

# Directory for storing shared data, like ADAS reaction rates and radiation fits
GKYL_SHARE_DIR ?= "${INSTALL_PREFIX}/gkylzero/share"
CFLAGS += -DGKYL_SHARE_DIR=$(GKYL_SHARE_DIR)

# Read MPI paths and flags if needed 
USING_MPI =
MPI_RPATH = 
MPI_INC_DIR = zero # dummy
MPI_LIB_DIR = .
ifeq (${USE_MPI}, 1)
	USING_MPI = yes
	MPI_INC_DIR = ${CONF_MPI_INC_DIR}
	MPI_LIB_DIR = ${CONF_MPI_LIB_DIR}
ifdef USING_NVCC
	MPI_RPATH = -Xlinker "-rpath,${CONF_MPI_LIB_DIR}"
else
	MPI_RPATH = -Wl,-rpath,${CONF_MPI_LIB_DIR}
endif
	MPI_LIBS = -lmpi
	CFLAGS += -DGKYL_HAVE_MPI
endif

# Read NCCL paths and flags if needed (needs MPI and NVCC)
USING_NCCL =
NCCL_INC_DIR = zero # dummy
NCCL_LIB_DIR = .
ifeq (${USE_NCCL}, 1)
ifdef USING_MPI
ifdef USING_NVCC
	USING_NCCL = yes
	NCCL_INC_DIR = ${CONF_NCCL_INC_DIR}
	NCCL_LIB_DIR = ${CONF_NCCL_LIB_DIR}
	NCCL_LIBS = -lnccl
	CFLAGS += -DGKYL_HAVE_NCCL
endif
endif
endif

# Read CUDSS paths and flags if needed (needs MPI and NVCC)
USING_CUDSS =
CUDSS_INC_DIR = zero # dummy
CUDSS_LIB_DIR = .
CUDSS_RPATH =
ifeq (${USE_CUDSS}, 1)
ifdef USING_NVCC
	USING_CUDSS = yes
	CUDSS_INC_DIR = ${CONF_CUDSS_INC_DIR}
	CUDSS_LIB_DIR = ${CONF_CUDSS_LIB_DIR}
ifdef USING_NVCC
	CUDSS_RPATH = -Xlinker "-rpath,${CONF_CUDSS_LIB_DIR}"
else
	CUDSS_RPATH = -Wl,-rpath,${CONF_CUDSS_LIB_DIR}
endif
	CUDSS_LIBS = -lcudss
	CFLAGS += -DGKYL_HAVE_CUDSS
endif
endif

# Read LUA paths and flags if needed 
USING_LUA =
LUA_RPATH = 
LUA_INC_DIR = zero # dummy
LUA_LIB_DIR = .
ifeq (${USE_LUA}, 1)
	USING_LUA = yes
	LUA_INC_DIR = ${CONF_LUA_INC_DIR}
	LUA_LIB_DIR = ${CONF_LUA_LIB_DIR}
	LUA_RPATH = -Wl,-rpath,${CONF_LUA_LIB_DIR}
	LUA_LIBS = -l${CONF_LUA_LIB}
	CFLAGS += -DGKYL_HAVE_LUA
endif

# Read ADAS paths and flags if needed 
USING_ADAS =
ADAS_INC_DIR = zero # dummy
ADAS_LIB_DIR = .
ifeq (${USE_ADAS}, 1)
	USING_ADAS = yes
	CFLAGS += -DGKYL_HAVE_ADAS
endif

# Build directory
ifdef USING_NVCC
	BUILD_DIR ?= cuda-build
else	
	BUILD_DIR ?= build
endif

# On OSX we should use Accelerate framework
ifeq ($(UNAME), Darwin)
	LAPACK_LIB_DIR = .
	LAPACK_INC = zero # dummy
	LAPACK_LIB = -framework Accelerate
	CFLAGS += -DGKYL_USING_FRAMEWORK_ACCELERATE
	SHFLAGS += -dynamiclib 
else
	SHFLAGS += -shared
endif

# For install shared-lib we need to pass extra flag to Mac
# compiler. See note below for ZERO_SH_INSTALL_LIB target.
SHFLAGS_INSTALL = ${SHFLAGS}
ifeq ($(UNAME), Darwin)
	SHFLAGS_INSTALL = ${SHFLAGS} -install_name ${PREFIX}/gkylzero/lib/libgkylzero.so
endif

# Header files 
HEADERS := $(wildcard minus/*.h) $(wildcard zero/*.h) $(wildcard apps/*.h) $(wildcard amr/*.h) $(wildcard kernels/*/*.h) $(wildcard data/adas/*.h)
# Headers to install
INSTALL_HEADERS := $(shell ls apps/gkyl_*.h zero/gkyl_*.h  amr/gkyl_*.h | grep -v "priv" | sort)
INSTALL_HEADERS += $(shell ls minus/*.h)

# all includes
INCLUDES = -Iminus -Iminus/STC/include -Izero -Iapps -Iamr -Iregression -I${BUILD_DIR} ${KERN_INCLUDES} -I${LAPACK_INC} -I${SUPERLU_INC} -I${MPI_INC_DIR} -I${NCCL_INC_DIR} -I${CUDSS_INC_DIR} -I${LUA_INC_DIR}

# Directories containing source code
SRC_DIRS := minus zero apps amr kernels data/adas

# List of regression and unit test
REGS := $(patsubst %.c,${BUILD_DIR}/%,$(wildcard regression/rt_*.c))
UNITS := $(patsubst %.c,${BUILD_DIR}/%,$(wildcard unit/ctest_*.c))
MPI_UNITS := $(patsubst %.c,${BUILD_DIR}/%,$(wildcard unit/mctest_*.c))
LUA_UNITS := $(patsubst %.c,${BUILD_DIR}/%,$(wildcard unit/lctest_*.c))
CI := $(patsubst %.c,${BUILD_DIR}/%,$(wildcard ci/*.c))

# list of includes from kernels
KERN_INC_DIRS = $(shell find $(SRC_DIRS) -type d)
KERN_INCLUDES = $(addprefix -I,$(KERN_INC_DIRS))

# We need to build CUDA unit-test objects
UNIT_CU_SRCS =
UNIT_CU_OBJS =
# There is some problem with the Vlasov and Maxwell kernels that is causing some unit builds to fail
ifdef USING_NVCC
#	UNIT_CU_SRCS = $(shell find unit -name *.cu)
	UNIT_CU_SRCS = unit/ctest_cusolver.cu unit/ctest_alloc_cu.cu unit/ctest_basis_cu.cu unit/ctest_array_cu.cu unit/ctest_mom_vlasov_cu.cu unit/ctest_range_cu.cu unit/ctest_rect_grid_cu.cu unit/ctest_wave_geom_cu.cu unit/ctest_wv_euler_cu.cu
ifdef USING_CUDSS
	UNIT_CU_SRCS += unit/ctest_cudss.cu
endif
	UNIT_CU_OBJS = $(UNIT_CU_SRCS:%=$(BUILD_DIR)/%.o)
endif

# List of link directories and libraries for unit and regression tests
EXEC_LIB_DIRS = -L${SUPERLU_LIB_DIR} -L${LAPACK_LIB_DIR} -L${BUILD_DIR} -L${MPI_LIB_DIR} -L${NCCL_LIB_DIR} -L${CUDSS_LIB_DIR} -L${LUA_LIB_DIR}
EXEC_EXT_LIBS = -lsuperlu ${LAPACK_LIB} ${CUDA_LIBS} ${MPI_RPATH} ${MPI_LIBS} ${NCCL_LIBS} ${CUDSS_RPATH} ${CUDSS_LIBS} ${LUA_RPATH} ${LUA_LIBS} -lm -lpthread -ldl
EXEC_LIBS = ${BUILD_DIR}/libgkylzero.so ${EXEC_EXT_LIBS}
EXEC_INSTALLED_LIBS = ${G0_RPATH} -lgkylzero ${EXEC_EXT_LIBS}
EXEC_RPATH = 

# Rpath for use in glua exectuable
G0_LIB_DIR = ${INSTALL_PREFIX}/gkylzero/lib
G0_RPATH = -Wl,-rpath,${G0_LIB_DIR}

# Build commands for C source
$(BUILD_DIR)/%.c.o: %.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(ARCH_FLAGS) $(INCLUDES) -c $< -o $@

# Build commands for CUDA source
$(BUILD_DIR)/%.cu.o: %.cu
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

# Unit tests
${BUILD_DIR}/unit/%: unit/%.c ${BUILD_DIR}/libgkylzero.so ${UNIT_CU_OBJS}
	$(MKDIR_P) ${BUILD_DIR}/unit
	${CC} ${CFLAGS} ${LDFLAGS} -o $@ $< -I. $(INCLUDES) ${UNIT_CU_OBJS} ${EXEC_LIB_DIRS} ${EXEC_RPATH} ${EXEC_LIBS}

# Regression tests
${BUILD_DIR}/regression/%: regression/%.c ${BUILD_DIR}/libgkylzero.so
	$(MKDIR_P) ${BUILD_DIR}/regression
	${CC} ${CFLAGS} ${LDFLAGS} -o $@ $< -I. $(INCLUDES) ${EXEC_LIB_DIRS} ${EXEC_RPATH} ${EXEC_LIBS}

# Automated regression system
${BUILD_DIR}/ci/%: ci/%.c ${BUILD_DIR}/libgkylzero.so
	$(MKDIR_P) ${BUILD_DIR}/ci
	${CC} ${CFLAGS} ${LDFLAGS} -o $@ $< -I. $(INCLUDES) ${EXEC_LIB_DIRS} ${EXEC_RPATH} ${EXEC_LIBS}

# Lua interpreter for testing Lua regression tests
${BUILD_DIR}/glua: regression/glua.c ${BUILD_DIR}/libgkylzero.so
	$(MKDIR_P) ${BUILD_DIR}
	${CC} ${CFLAGS} ${LDFLAGS} -o $@ $< -I. $(INCLUDES) ${EXEC_LIB_DIRS} ${EXEC_RPATH} ${EXEC_LIBS}

# Lua interpreter for testing Lua regression tests
${BUILD_DIR}/glua-install: regression/glua.c ${BUILD_DIR}/libgkylzero.so
	$(MKDIR_P) ${BUILD_DIR}
	${CC} ${CFLAGS} ${LDFLAGS} -o $@ $< -I. $(INCLUDES) ${EXEC_LIB_DIRS} ${EXEC_INSTALLED_LIBS} 

# Amalgamated header file
${BUILD_DIR}/gkylzero.h:
	$(MKDIR_P) ${BUILD_DIR}
	./minus/gengkylzeroh.sh > ${BUILD_DIR}/gkylzero.h

# Specialized build commands for kernels when using nvcc
ifdef USING_NVCC

# Unfortunately, due to the limitations of the NVCC compiler to treat
# device code in C files, we need to force compile the kernel code
# using the -x cu flag

$(BUILD_DIR)/kernels/advection/%.c.o : kernels/advection/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/bgk/%.c.o : kernels/bgk/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/bin_op/%.c.o : kernels/bin_op/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/canonical_pb/%.c.o : kernels/canonical_pb/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/dg_diffusion_fluid/%.c.o : kernels/dg_diffusion_fluid/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/dg_diffusion_vlasov/%.c.o : kernels/dg_diffusion_vlasov/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/dg_diffusion_gyrokinetic/%.c.o : kernels/dg_diffusion_gyrokinetic/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/dg_diffusion_gen/%.c.o : kernels/dg_diffusion_gen/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/euler/%.c.o : kernels/euler/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/gyrokinetic/%.c.o : kernels/gyrokinetic/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/lbo/%.c.o : kernels/lbo/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/maxwell/%.c.o : kernels/maxwell/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/pkpm/%.c.o : kernels/pkpm/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/prim_vars/%.c.o : kernels/prim_vars/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/rad/%.c.o : kernels/rad/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/vlasov/%.c.o : kernels/vlasov/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/sr_vlasov/%.c.o : kernels/sr_vlasov/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/basis/%.c.o : kernels/basis/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/fem_poisson/%.c.o : kernels/fem_poisson/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/fem_parproj/%.c.o : kernels/fem_parproj/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/fem_poisson_perp/%.c.o : kernels/fem_poisson_perp/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/ambi_bolt_potential/%.c.o : kernels/ambi_bolt_potential/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/array_integrate/%.c.o : kernels/array_integrate/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/deflate_zsurf/%.c.o : kernels/deflate_zsurf/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/positivity_shift/%.c.o : kernels/positivity_shift/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

endif

## GkylZero Library 
ZERO := libgkylzero
SRCS := $(shell find $(SRC_DIRS) -name *.c)
ifdef USING_NVCC
	SRCS += $(shell find zero -name *.cu)
endif
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

ZERO_SH_LIB := $(BUILD_DIR)/$(ZERO).so
$(ZERO_SH_LIB): $(OBJS)
	$(MKDIR_P) $(dir $@)
	${CC} ${SHFLAGS} ${LDFLAGS} ${OBJS} ${EXEC_LIB_DIRS} ${EXEC_EXT_LIBS} -o $@

# Due to an issue with shared-lib linking on the Mac, we need to build
# a separate shared lib to install. This one has the install path
# hard-coded into the library itself, so external execs, like gkyl,
# can link to the library properly. Perhaps there is a another way to
# do this, don't know. -- AH, Feb 4th 2023.
ZERO_SH_INSTALL_LIB := $(BUILD_DIR)/$(ZERO)-install.so
$(ZERO_SH_INSTALL_LIB): $(OBJS)
	$(MKDIR_P) $(dir $@)
	${CC} ${SHFLAGS_INSTALL} ${LDFLAGS} ${OBJS} ${EXEC_LIB_DIRS} ${EXEC_EXT_LIBS} -o $@

## All libraries build targets completed at this point

.PHONY: all
all: ${BUILD_DIR}/gkylzero.h ${ZERO_SH_LIB} ## Build libraries and amalgamated header
	${MKDIR_P} ${INSTALL_PREFIX}/gkylzero/share/adas
	cp ./data/adas/radiation_fit_parameters.txt ${INSTALL_PREFIX}/gkylzero/share/adas

# Explicit targets to build unit and regression tests
unit: ${ZERO_SH_LIB} ${UNITS} ${MPI_UNITS} ${LUA_UNITS} ## Build unit tests
regression: ${ZERO_SH_LIB} ${REGS} regression/rt_arg_parse.h ${BUILD_DIR}/glua ## Build regression tests
ci: ${ZERO_SH_LIB} ${CI} ## Build automated regression system
glua: ${BUILD_DIR}/glua ## Build Lua interpreter

.PHONY: check mpicheck
# Run all unit tests
check: ${UNITS} ## Build (if needed) and run all unit tests
	$(foreach unit,${UNITS},echo $(unit); $(unit) -E;)

# Run all unit tests needing MPI
mpicheck: ${MPI_UNITS} ## Build (if needed) and run all unit tests needing MPI
	$(foreach unit,${MPI_UNITS},echo $(unit); $(unit) -E -M;)

# Set full paths to install location so config.mak is used from there
G0_SHARE_INSTALL_PREFIX=${INSTALL_PREFIX}/gkylzero/share
SED_REPS_STR=s,G0_SHARE_INSTALL_PREFIX_TAG,${G0_SHARE_INSTALL_PREFIX},g

install: all $(ZERO_SH_INSTALL_LIB) ${BUILD_DIR}/glua ## Install library and headers
# Construct install directories
	$(MKDIR_P) ${INSTALL_PREFIX}/gkylzero/include
	${MKDIR_P} ${INSTALL_PREFIX}/gkylzero/lib
	${MKDIR_P} ${INSTALL_PREFIX}/gkylzero/bin
	${MKDIR_P} ${INSTALL_PREFIX}/gkylzero/share
	${MKDIR_P} ${INSTALL_PREFIX}/gkylzero/share/adas
# Headers
	cp ${INSTALL_HEADERS} ${INSTALL_PREFIX}/gkylzero/include
	./minus/gengkylzeroh.sh > ${INSTALL_PREFIX}/gkylzero/include/gkylzero.h
# libraries
#	strip ${ZERO_SH_INSTALL_LIB}  ## MF 2023/10/26: Causes a problem in Macs.
	cp -f ${ZERO_SH_INSTALL_LIB} ${INSTALL_PREFIX}/gkylzero/lib/libgkylzero.so
# Examples
	test -e config.mak && cp -f config.mak ${INSTALL_PREFIX}/gkylzero/share/config.mak || echo "No config.mak"
	sed ${SED_REPS_STR} Makefile.sample > ${INSTALL_PREFIX}/gkylzero/share/Makefile
	cp -f regression/rt_arg_parse.h ${INSTALL_PREFIX}/gkylzero/include/rt_arg_parse.h
	cp -f regression/rt_vlasov_twostream_p2.c ${INSTALL_PREFIX}/gkylzero/share/rt_vlasov_twostream_p2.c

install-glua: install ${BUILD_DIR}/glua-install ## Install the glua exectuable
# glua executable
	cp -f ${BUILD_DIR}/glua-install ${INSTALL_PREFIX}/gkylzero/bin/glua


.PHONY: clean
clean: ## Clean build output
	rm -rf ${BUILD_DIR}

.PHONY: cleanur
cleanur: ## Delete the unit and regression test executables
	rm -rf ${BUILD_DIR}/unit ${BUILD_DIR}/regression ${BUILD_DIR}/ci ${BUILD_DIR}/glua

.PHONY: cleanr
cleanr: ## Delete the regression test executables
	rm -rf ${BUILD_DIR}/regression

.PHONY: cleanu
cleanu: ## Delete the unit test executables
	rm -rf ${BUILD_DIR}/unit

.PHONY: cleanc
cleanc: ## Delete the automated regression test executables
	rm -rf ${BUILD_DIR}/ci

# include dependencies
-include $(DEPS)

# command to make dir
MKDIR_P ?= mkdir -p

# From: https://www.client9.com/self-documenting-makefiles/
.PHONY: help
help: ## Show help
	@echo "GkylZero Makefile help. You can set parameters on the command line:"
	@echo ""
	@echo "make CC=nvcc -j"
	@echo ""
	@echo "Or run the configure script to set various parameters. Usually"
	@echo "defaults are all you need, specially if the dependencies are in"
	@echo "${HOME}/gkylsoft and you are using standard compilers (not building on GPUs)"
	@echo ""
	@echo "See ./configure --help for usage of configure script"
	@echo ""
	@awk -F ':|##' '/^[^\t].+?:.*?##/ {\
        printf "\033[36m%-30s\033[0m %s\n", $$1, $$NF \
        }' $(MAKEFILE_LIST)
