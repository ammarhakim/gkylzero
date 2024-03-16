#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h> 


struct gkyl_mj_moments
{
  struct gkyl_rect_grid grid;
  struct gkyl_basis conf_basis, phase_basis;

  // struct gkyl_mom_calc *ncalc;
  struct gkyl_dg_updater_moment *ncalc; 
  struct gkyl_dg_updater_moment *vbicalc;
  struct gkyl_dg_updater_moment *Pcalc;
  struct gkyl_array *num_ratio; 
  struct gkyl_array *num_vb;  
  struct gkyl_array *vb_dot_nvb;  
  struct gkyl_array *n_minus_vb_dot_nvb;  
  struct gkyl_array *V_drift;   
  struct gkyl_array *gamma;
  struct gkyl_array *GammaV2;
  struct gkyl_array *Gamma_inv;
  struct gkyl_array *pressure;
  struct gkyl_array *temperature;

  struct gkyl_dg_bin_op_mem *mem;
};