#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

struct gkyl_correct_vlasov_lte
{
  struct gkyl_rect_grid grid;
  struct gkyl_basis conf_basis, phase_basis;

  struct gkyl_array *num_ratio;  
  struct gkyl_array *num_vb;   
  struct gkyl_array *V_drift;  

  struct gkyl_dg_bin_op_mem *mem;     
  struct gkyl_array *n_stationary, *vbi, *T_stationary;
  struct gkyl_array *dn, *dvbi, *dT;
  struct gkyl_array *ddn, *ddvbi, *ddT;

  struct gkyl_maxwellian_moments *moments_up;
  struct gkyl_array *moms;
  struct gkyl_array *moms_target;
  struct gkyl_proj_mj_on_basis *proj_mj;
  struct gkyl_proj_maxwellian_on_basis *proj_max;

  enum gkyl_model_id model_id;

  // error estimate n, vb, T, 0 - success., num. picard iterations
  double error_n; 
  double error_vb[3];
  double error_T;
  int status; 
  int niter;
};