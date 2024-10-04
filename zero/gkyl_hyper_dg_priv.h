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
