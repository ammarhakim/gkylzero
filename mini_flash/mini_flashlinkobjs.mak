# -*- makefile-gmake -*-

# mini_flash link objects

mini_flash_LINK_OBJS = ../$(BUILD_DIR)/mini_flash/zero/*.o ../$(BUILD_DIR)/mini_flash/apps/*.o ../$(BUILD_DIR)/mini_flash/$(KERNELS_DIR)/*/*.o
# ... with GPUs
mini_flash_LINK_CU_OBJS = 
