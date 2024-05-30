// Private header: not for direct use
#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

struct gkyl_vlasov_lte_correct
{
  struct gkyl_basis conf_basis; // Configuration-space basis
  struct gkyl_basis phase_basis; // Phase-space basis
  int num_conf_basis; // Number of configuration-space basis functions
  int vdim; // Number of velocity dimensions
  enum gkyl_model_id model_id; // Enum identifier for model type (e.g., SR, see gkyl_eqn_type.h) 

  struct gkyl_array *moms_iter;
  struct gkyl_array *d_moms;
  struct gkyl_array *dd_moms;

  struct gkyl_vlasov_lte_moments *moments_up;
  struct gkyl_vlasov_lte_proj_on_basis *proj_lte;

  // error estimate (n, V_drift, T/m), 0 - success., num. picard iterations
  double error[5]; // absolute value of difference in cell averages between iteration and target
  double eps; // tolerance for the iterator
  int max_iter; // number of total iterations
  
  bool use_gpu; // Boolean if we are performing projection on device.
  double *error_cu; // error on device if using GPUs 
  struct gkyl_array *abs_diff_moms;

  struct gkyl_velocity_map *vel_map; // Velocity space mapping object.
};
