# -*- makefile-gmake -*-

# Gyrokinetic link objects

GYROKINETIC_LINK_OBJS = ../$(BUILD_DIR)/gyrokinetic/zero/*.o ../$(BUILD_DIR)/gyrokinetic/apps/*.o ../$(BUILD_DIR)/gyrokinetic/data/adas/*.o ../$(BUILD_DIR)/gyrokinetic/kernels/*/*.o
# ... with GPUs
GYROKINETIC_LINK_CU_OBJS = 