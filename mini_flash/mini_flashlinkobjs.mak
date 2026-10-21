# -*- makefile-gmake -*-

# mini_flash link objects

MINI_FLASH_LINK_OBJS = ../$(BUILD_DIR)/mini_flash/zero/*.o ../$(BUILD_DIR)/mini_flash/apps/*.o ../$(BUILD_DIR)/mini_flash/$(KERNELS_DIR)/*/*.o
# ... with GPUs
MINI_FLASH_LINK_CU_OBJS = 
