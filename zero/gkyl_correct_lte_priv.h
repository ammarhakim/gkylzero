#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

struct gkyl_correct_vlasov_lte
{
  struct gkyl_basis conf_basis; // Configuration-space basis
  struct gkyl_basis phase_basis; // Phase-space basis
  int num_conf_basis; // Number of configuration-space basis functions
  int vdim; // Number of velocity dimensions
  enum gkyl_model_id model_id; // Enum identifier for model type (e.g., SR, see gkyl_eqn_type.h)

  struct gkyl_array *moms_dens_corr;
  struct gkyl_array *num_ratio; 
  struct gkyl_dg_bin_op_mem *mem;  

  struct gkyl_array *moms_iter;
  struct gkyl_array *d_moms;
  struct gkyl_array *dd_moms;

  struct gkyl_maxwellian_moments *moments_up;
  struct gkyl_proj_mj_on_basis *proj_mj;
  struct gkyl_proj_maxwellian_on_basis *proj_max;

  // error estimate n, vb, T, 0 - success., num. picard iterations
  double error_n; 
  double error_vb[3];
  double error_T;
};
