#pragma once

#include <gkyl_dg_eqn.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

struct gkyl_hyper_dg {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  int num_basis; // number of basis functions
  int num_up_dirs; // number of update directions
  int update_dirs[GKYL_MAX_DIM]; // directions to update
  // zero_flux_flags[d] == 1 means zero-flux BC in 'd'
  int zero_flux_flags[GKYL_MAX_DIM];
  int update_vol_term; // should we update volume term?
  const struct gkyl_dg_eqn *equation; // equation object

  uint32_t flags;
  struct gkyl_hyper_dg *on_dev; // pointer to itself or device data
};

GKYL_CU_DH
static void
create_offsets(const int num_up_dirs, const int update_dirs[2], 
  const struct gkyl_range *range, const int idxc[GKYL_MAX_DIM], long offsets[9])
{
  
  // Check if we're at an upper or lower edge in each direction
  bool is_edge_upper[2], is_edge_lower[2];
  for (int i=0; i<num_up_dirs; ++i) {
    is_edge_lower[i] = idxc[update_dirs[i]] == range->lower[update_dirs[i]];
    is_edge_upper[i] = idxc[update_dirs[i]] == range->upper[update_dirs[i]];
  }

  // Construct the offsets *only* in the directions being updated.
  // No need to load the neighbors that are not needed for the update.
  int lower_offset[GKYL_MAX_DIM] = {0};
  int upper_offset[GKYL_MAX_DIM] = {0};
  for (int d=0; d<num_up_dirs; ++d) {
    int dir = update_dirs[d];
    lower_offset[dir] = -1 + is_edge_lower[d];
    upper_offset[dir] = 1 - is_edge_upper[d];
  }  

  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, lower_offset, upper_offset);
  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);
  // construct list of offsets
  int count = 0;
  while (gkyl_range_iter_next(&iter3))
    offsets[count++] = gkyl_range_offset(range, iter3.idx);

}

GKYL_CU_DH
static int 
idx_to_inloup_ker(int dim, const int *idx, const int *dirs, const int *num_cells) {
  int iout = 0;

  for (int d=0; d<dim; ++d) {
    if (idx[dirs[d]] == 1) {
      iout = 2*iout+(int)(pow(3,d)+0.5);
    } else if (idx[dirs[d]] == num_cells[dirs[d]]) {
      iout = 2*iout+(int)(pow(3,d)+0.5)+1;
    }
  }
  return iout;
}
