# -*- makefile-gmake -*-

# Core link objects

CORE_LINK_OBJS = ../$(BUILD_DIR)/core/minus/*.o ../$(BUILD_DIR)/core/kernels/*/*.o ../$(BUILD_DIR)/core/zero/*.o ../$(BUILD_DIR)/core/lua/*/*.o ../$(BUILD_DIR)/core/apps/*.o
# ... with GPUs
CORE_LINK_CU_OBJS = ../$(BUILD_DIR)/core/unit/*.cu.o
