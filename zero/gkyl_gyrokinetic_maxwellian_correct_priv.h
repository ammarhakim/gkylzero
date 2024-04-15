// Private header: not for direct use
#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

struct gkyl_gyrokinetic_maxwellian_correct
{
  struct gkyl_basis conf_basis; // Configuration-space basis
  struct gkyl_basis phase_basis; // Phase-space basis
  int num_conf_basis; // Number of configuration-space basis functions
  const struct gk_geometry *gk_geom; // Geometry struct
  bool divide_jacobgeo; // Boolean for if we are dividing out the configuration-space Jacobian from density
  bool use_last_converged; // Boolean for if we are using the results of the iterative scheme
                           // *even if* the scheme fails to converge. 
  double mass; // Species mass

  struct gkyl_array *moms_iter;
  struct gkyl_array *d_moms;
  struct gkyl_array *dd_moms;
  struct gkyl_array *num_ratio; // num_ratio = n/n0 for correcting density
  struct gkyl_dg_bin_op_mem *mem; // bin_op memory for correcting density

  struct gkyl_gyrokinetic_maxwellian_moments *moments_up;
  struct gkyl_proj_maxwellian_on_basis *proj_max_prim;

  // error estimate (n, u_par, T/m), 0 - success., num. picard iterations
  double error[3]; // absolute value of difference in cell averages between iteration and target
  double eps; // tolerance for the iterator
  int max_iter; // number of total iterations
  
  bool use_gpu; // Boolean if we are performing projection on device.
  double *error_cu; // error on device if using GPUs 
  struct gkyl_array *abs_diff_moms;
};
