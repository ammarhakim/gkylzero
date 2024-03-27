// Private header: not for direct use
#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>

struct gkyl_vlasov_lte_moments
{
  struct gkyl_basis conf_basis; // Configuration-space basis
  struct gkyl_basis phase_basis; // Phase-space basis
  int num_conf_basis; // Number of configuration-space basis functions
  int vdim; // Number of velocity dimensions
  enum gkyl_model_id model_id; // Enum identifier for model type (e.g., SR, see gkyl_eqn_type.h)
  double mass; // Species mass

  struct gkyl_array *M0; 
  struct gkyl_array *M1i;  
  struct gkyl_array *V_drift;
  struct gkyl_array *V_drift_dot_M1i;  
  struct gkyl_array *pressure;
  struct gkyl_array *temperature;
  struct gkyl_dg_bin_op_mem *mem;

  union {
    // special relativistic Vlasov-Maxwell model
    struct {
      struct gkyl_array *GammaV2;
      struct gkyl_array *GammaV_inv;
      struct gkyl_array *M0_minus_V_drift_dot_M1i;  
    };
  };

  struct gkyl_dg_updater_moment *M0_calc; 
  struct gkyl_dg_updater_moment *M1i_calc;
  struct gkyl_dg_updater_moment *Pcalc;
};
