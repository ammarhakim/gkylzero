# -*- makefile-gmake -*-

# Vlasov link objects

VLASOV_LINK_OBJS = ../$(BUILD_DIR)/vlasov/kernels/*/*.o ../$(BUILD_DIR)/vlasov/zero/*.o ../$(BUILD_DIR)/vlasov/apps/*.o
# ... with GPUs
VLASOV_LINK_CU_OBJS = 