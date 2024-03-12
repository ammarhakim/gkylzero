#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_correct_mj gkyl_correct_mj;

struct gkyl_correct_mj
{
  struct gkyl_rect_grid grid;
  struct gkyl_basis conf_basis, phase_basis;

  struct gkyl_dg_updater_moment *m0calc;  
  struct gkyl_dg_updater_moment *m1icalc; 
  struct gkyl_array *num_ratio;  
  struct gkyl_array *num_vb;   
  struct gkyl_array *V_drift;  
  struct gkyl_array *gamma;   
  struct gkyl_array *vb_dot_nvb; 
  struct gkyl_array *n_minus_vb_dot_nvb; 

  struct gkyl_dg_bin_op_mem *mem;     
  struct gkyl_array *m0, *m1i, *m2;
  struct gkyl_array *dm0, *dm1i, *dm2;
  struct gkyl_array *ddm0, *ddm1i, *ddm2;

  struct gkyl_mj_moments *mj_moms;
  struct gkyl_proj_mj_on_basis *proj_mj;

  // error estimate n, vb, T, 0 - success., num. picard iterations
  double error_n; 
  double error_vb[3];
  double error_T;
  int status; 
  int niter;
};

/**
 * Create new updater to correct a Maxwellian to match specified
 * moments.
 *
 * @param grid Grid on which updater lives
 * @param conf_basis Conf space basis functions
 * @param phase_basis Phase space basis functions
 * @param conf_range configuration space range
 * @param vel_range velocity space range
 * @param conf_local_cells Number of cells in local config-space
 * @param conf_local_ext_cells Number of cells in local extended config-space
 */
gkyl_correct_mj *gkyl_correct_mj_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *conf_range_ext, const struct gkyl_range *vel_range, 
  const struct gkyl_array *p_over_gamma, const struct gkyl_array *gamma, const struct gkyl_array *gamma_inv, 
  bool use_gpu);

/**
 * Fix the Maxwell-Juttner so that it's moments match desired m0 moment.
 *
 * @param cmj MJ correction updater
 * @param p_over_gamma velocity array
 * @param fout Distribution function to fix (modified in-place)
 * @param m0 Desired lab-frame number density
 * @param m1i specified velocity of f in the stationary frame
 * @param phase_local Local phase-space range
 * @param conf_local Local configuration space range
 */
void gkyl_correct_mj_fix_m0(gkyl_correct_mj *cmj, 
  struct gkyl_array *fout,
  const struct gkyl_array *m0, const struct gkyl_array *m1i,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local);

/**
 * Fix the Maxwell-Juttner so that it's moments match desired moments.
 * NOTE: If this algorithm fails, the returns the original distribution function
 * with only the m0 moment corrected (i.e. runs: gkyl_correct_mj_fix_m0())
 *
 * @param cmj MJ correction updater
 * @param fout Distribution function to fix (modified in-place)
 * @param m0 Desired lab-frame number density
 * @param m1i Desired lab-frame velocity
 * @param m2 Desired lab-frame temperature
 * @param phase_local Local phase-space range
 * @param conf_local Local configuration space range
 * @param vel_basis Vel space basis functions
 * @param vel_grid Vel-Grid on which updater lives
 */
void gkyl_correct_mj_fix(gkyl_correct_mj *cmj,
  struct gkyl_array *fout,
  const struct gkyl_array *m0, const struct gkyl_array *m1i, const struct gkyl_array *m2,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local,
  int poly_order);

/**
 * Delete updater.
 *
 * @param cmj Updater to delete.
 */
void gkyl_correct_mj_release(gkyl_correct_mj *cmj);
