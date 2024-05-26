#pragma once

#include <gkyl_bc_excision.h>

// Private header for bc_excision, not for direct use in user code.

// Primary struct in this updater.
struct gkyl_bc_excision {
  int tan_dir; // Tangential direction to transport the field in.
  int tan_dir_num_cellsD2; // Number of cells in tangential direction / 2.
  int num_basis; // Number of basis representing the field.
  struct gkyl_range ghost_r; // Ghost range.
  struct gkyl_range buff_r; // Buffer range.
  bool use_gpu; // Whether to apply the BC with the GPU.
};

#ifdef GKYL_HAVE_CUDA

/**
 * Apply the excision BC on the GPU.
 *
 * @param up Excision BC updater.
 * @param ghost_buffer Buffer with the values of the corresponding ghost cells.
 * @param distf Distribution function to apply BC to.
 */
void
gkyl_bc_excision_advance_cu(const struct gkyl_bc_excision *up,
  const struct gkyl_array *ghost_buffer, struct gkyl_array *distf);

#endif
