# -*- makefile-gmake -*-

# Type make help to see help for this Makefile

ARCH_FLAGS ?= -march=native
CUDA_ARCH ?= 70
# Warning flags: -Wall -Wno-unused-variable -Wno-unused-function -Wno-missing-braces
CFLAGS ?= -O3 -g -ffast-math -fPIC -MMD -MP
LDFLAGS = 
PREFIX ?= ${HOME}/gkylsoft

# determine OS we are running on
UNAME = $(shell uname)

# Default lapack include and libraries: we prefer linking to static library
LAPACK_INC = ${HOME}/gkylsoft/OpenBLAS/include
LAPACK_LIB_DIR = ${HOME}/gkylsoft/OpenBLAS/lib
LAPACK_LIB = -lopenblas

# SuperLU includes and librararies
SUPERLU_INC = ${HOME}/gkylsoft/superlu/include
ifeq ($(UNAME_S),Linux)
	SUPERLU_LIB_DIR = ${HOME}/gkylsoft/superlu/lib64
	SUPERLU_LIB = ${HOME}/gkylsoft/superlu/lib64/libsuperlu.a
else
	SUPERLU_LIB_DIR = ${HOME}/gkylsoft/superlu/lib
	SUPERLU_LIB = ${HOME}/gkylsoft/superlu/lib/libsuperlu.a
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

# Header files
HEADERS := $(wildcard minus/*.h) $(wildcard zero/*.h) $(wildcard apps/*.h) $(wildcard kernels/*/*.h)
# Headers to install
INSTALL_HEADERS := $(shell ls apps/gkyl_*.h zero/gkyl_*.h | grep -v "priv" | sort)
INSTALL_HEADERS += $(shell ls minus/*.h)

# all includes
INCLUDES = -Iminus -Iminus/STC/include -Izero -Iapps -Iregression -I${BUILD_DIR} ${KERN_INCLUDES} -I${LAPACK_INC} -I${SUPERLU_INC}

# Directories containing source code
SRC_DIRS := minus zero apps kernels

# List of regression and unit tests
REGS := $(patsubst %.c,${BUILD_DIR}/%,$(wildcard regression/rt_*.c))
UNITS := $(patsubst %.c,${BUILD_DIR}/%,$(wildcard unit/ctest_*.c))

# list of includes from kernels
KERN_INC_DIRS = $(shell find $(SRC_DIRS) -type d)
KERN_INCLUDES = $(addprefix -I,$(KERN_INC_DIRS))

# List of link directories and libraries for unit and regression tests
EXEC_LIB_DIRS = -L${SUPERLU_LIB_DIR} -L${LAPACK_LIB_DIR} -L${BUILD_DIR}
EXEC_EXT_LIBS = -lsuperlu ${LAPACK_LIB} ${CUDA_LIBS} -lm -lpthread
EXEC_LIBS = -lgkylzero ${EXEC_EXT_LIBS}
EXEC_RPATH = -Wl,-rpath,${BUILD_DIR}

# Build commands for C source
$(BUILD_DIR)/%.c.o: %.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(ARCH_FLAGS) $(INCLUDES) -c $< -o $@

# Build commands for CUDA source
$(BUILD_DIR)/%.cu.o: %.cu
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

# Unit tests
${BUILD_DIR}/unit/%: unit/%.c ${ZERO_SH_LIB} ${UNIT_CU_OBJS} ${UNIT_CU_SRCS}
	$(MKDIR_P) ${BUILD_DIR}/unit
	${CC} ${CFLAGS} ${LDFLAGS} -o $@ $< -I. $(INCLUDES) ${EXEC_LIB_DIRS} ${EXEC_RPATH} ${EXEC_LIBS}

# Regression tests
${BUILD_DIR}/regression/%: regression/%.c ${ZERO_SH_LIB}
	$(MKDIR_P) ${BUILD_DIR}/regression
	${CC} ${CFLAGS} ${LDFLAGS} -o $@ $< -I. $(INCLUDES) ${EXEC_LIB_DIRS} ${EXEC_RPATH} ${EXEC_LIBS}

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

$(BUILD_DIR)/kernels/bin_op/%.c.o : kernels/bin_op/%.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(NVCC_FLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_DIR)/kernels/dg_diffusion/%.c.o : kernels/dg_diffusion/%.c
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

$(BUILD_DIR)/kernels/vlasov/%.c.o : kernels/vlasov/%.c
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

## All libraries build targets completed at this point

.PHONY: all
all: ${BUILD_DIR}/gkylzero.h ${ZERO_SH_LIB} ## Build libraries and amalgamated header

# Explicit targets to build unit and regression tests
unit: ${UNITS} ## Build unit tests
regression: ${REGS} ## Build regression tests

.PHONY: check
# Run all unit tests
check: ${UNITS} ## Build (if needed) and run all unit tests
	$(foreach unit,${UNITS},echo $(unit); $(unit) -E;)

install: all ## Install library and headers
# Construct install directories
	$(MKDIR_P) ${PREFIX}/gkylzero/include
	${MKDIR_P} ${PREFIX}/gkylzero/lib
	${MKDIR_P} ${PREFIX}/gkylzero/bin
	${MKDIR_P} ${PREFIX}/gkylzero/share
	${MKDIR_P} ${PREFIX}/gkylzero/scripts
# Headers
	cp ${INSTALL_HEADERS} ${PREFIX}/gkylzero/include
	./minus/gengkylzeroh.sh > ${PREFIX}/gkylzero/include/gkylzero.h
# libraries
	cp -f ${ZERO_SH_LIB} ${PREFIX}/gkylzero/lib
# Examples
	cp -f Makefile.sample ${PREFIX}/gkylzero/share/Makefile
	cp -f regression/rt_arg_parse.h ${PREFIX}/gkylzero/share/rt_arg_parse.h
	cp -f regression/rt_twostream.c ${PREFIX}/gkylzero/share/rt_twostream.c
# Lua wrappers
	cp -f inf/Vlasov.lua ${PREFIX}/gkylzero/lib/
	cp -f inf/Moments.lua ${PREFIX}/gkylzero/lib/
# Misc scripts
	cp -f scripts/*.sh ${PREFIX}/gkylzero/scripts

.PHONY: clean
clean: ## Clean build output
	rm -rf ${BUILD_DIR}

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
