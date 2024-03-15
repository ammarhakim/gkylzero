#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

struct gkyl_correct_mj
{
  struct gkyl_rect_grid grid;
  struct gkyl_basis conf_basis, phase_basis;

  struct gkyl_dg_updater_moment *ncalc;  
  struct gkyl_dg_updater_moment *vbicalc; 
  struct gkyl_array *num_ratio;  
  struct gkyl_array *num_vb;   
  struct gkyl_array *V_drift;  
  struct gkyl_array *gamma;   
  struct gkyl_array *vb_dot_nvb; 
  struct gkyl_array *n_minus_vb_dot_nvb; 

  struct gkyl_dg_bin_op_mem *mem;     
  struct gkyl_array *n, *vbi, *T;
  struct gkyl_array *dn, *dvbi, *dT;
  struct gkyl_array *ddn, *ddvbi, *ddT;

  struct gkyl_mj_moments *mj_moms;
  struct gkyl_proj_mj_on_basis *proj_mj;

  // error estimate n, vb, T, 0 - success., num. picard iterations
  double error_n; 
  double error_vb[3];
  double error_T;
  int status; 
  int niter;
};