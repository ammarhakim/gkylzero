# -*- makefile-gmake -*-

# Moments link objects

MOMENTS_LINK_OBJS = ../$(BUILD_DIR)/moments/zero/*.o ../$(BUILD_DIR)/moments/kernels/*/*.o ../$(BUILD_DIR)/moments/apps/*.o ../$(BUILD_DIR)/moments/amr/*.o
# ... with GPUs
MOMENTS_LINK_CU_OBJS = ../$(BUILD_DIR)/moments/unit/*.cu.o
