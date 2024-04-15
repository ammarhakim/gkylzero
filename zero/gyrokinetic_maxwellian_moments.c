#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_dg_updater_moment_gyrokinetic.h>
#include <gkyl_gyrokinetic_maxwellian_moments.h>
#include <gkyl_gyrokinetic_maxwellian_moments_priv.h>

struct gkyl_gyrokinetic_maxwellian_moments*
gkyl_gyrokinetic_maxwellian_moments_inew(const struct gkyl_gyrokinetic_maxwellian_moments_inp *inp)
{
  gkyl_gyrokinetic_maxwellian_moments *up = gkyl_malloc(sizeof(*up));

  up->conf_basis = *inp->conf_basis;
  up->phase_basis = *inp->phase_basis;
  up->num_conf_basis = inp->conf_basis->num_basis;
  // Determine factor to divide out of temperature computation
  // If 1x1v, up->vdim = 1, otherwise up->vdim = 3
  if (up->phase_basis.ndim - up->conf_basis.ndim == 1) {
    up->vdim = up->phase_basis.ndim - up->conf_basis.ndim;
  }
  else {
    up->vdim = 3;
  }
  up->gk_geom = gkyl_gk_geometry_acquire(inp->gk_geom);
  up->divide_jacobgeo = inp->divide_jacobgeo;

  long conf_range_ncells = inp->conf_range->volume;
  long conf_range_ext_ncells = inp->conf_range_ext->volume;

  if (inp->use_gpu) {
    up->M0 = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_range_ext_ncells);
    up->M1 = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_range_ext_ncells);
    up->u_par = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_range_ext_ncells);
    up->u_par_dot_M1 = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_range_ext_ncells);
    up->pressure = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_range_ext_ncells);
    up->temperature = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_range_ext_ncells);
    // Bin op memory for needed weak divisions
    up->mem = gkyl_dg_bin_op_mem_cu_dev_new(conf_range_ncells, up->num_conf_basis);
    if (up->vdim == 3) {
      // Additional moments if computing Bi-Maxwellian moments (Tpar, Tperp)
      up->p_par = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_range_ext_ncells);
      up->t_par = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_range_ext_ncells);    
      up->p_perp = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_range_ext_ncells);
      up->t_perp = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_range_ext_ncells); 
    }
  }
  else {
    up->M0 = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_range_ext_ncells);
    up->M1 = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_range_ext_ncells);
    up->u_par = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_range_ext_ncells);
    up->u_par_dot_M1 = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_range_ext_ncells);
    up->pressure = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_range_ext_ncells);
    up->temperature = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_range_ext_ncells);
    // Bin op memory for needed weak divisions
    up->mem = gkyl_dg_bin_op_mem_new(conf_range_ncells, up->num_conf_basis);
    if (up->vdim == 3) {
      // Additional moments if computing Bi-Maxwellian moments (Tpar, Tperp)
      up->p_par = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_range_ext_ncells);
      up->t_par = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_range_ext_ncells);
      up->p_perp = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_range_ext_ncells);
      up->t_perp = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_range_ext_ncells);
    }
  }

  // Moment calculator for needed moments (M0, M1, and M2)
  up->M0_calc = gkyl_dg_updater_moment_gyrokinetic_new(inp->phase_grid, inp->conf_basis,
    inp->phase_basis, inp->conf_range, inp->vel_range, inp->mass, inp->gk_geom, "M0", 0, inp->use_gpu);
  up->M1_calc = gkyl_dg_updater_moment_gyrokinetic_new(inp->phase_grid, inp->conf_basis,
    inp->phase_basis, inp->conf_range, inp->vel_range, inp->mass, inp->gk_geom, "M1", 0, inp->use_gpu);
  up->M2_calc = gkyl_dg_updater_moment_gyrokinetic_new(inp->phase_grid, inp->conf_basis,
    inp->phase_basis, inp->conf_range, inp->vel_range, inp->mass, inp->gk_geom, "M2", 0, inp->use_gpu);   

  if (up->vdim == 3) {
    // Additional moment calculators for Bi-Maxwellian moments (M2par, M2perp) 
    up->M2_par_calc = gkyl_dg_updater_moment_gyrokinetic_new(inp->phase_grid, inp->conf_basis,
      inp->phase_basis, inp->conf_range, inp->vel_range, inp->mass, inp->gk_geom, "M2par", 0, inp->use_gpu);   
    up->M2_perp_calc = gkyl_dg_updater_moment_gyrokinetic_new(inp->phase_grid, inp->conf_basis,
      inp->phase_basis, inp->conf_range, inp->vel_range, inp->mass, inp->gk_geom, "M2perp", 0, inp->use_gpu);  
  } 

  return up;
}

void 
gkyl_gyrokinetic_maxwellian_density_moment_advance(struct gkyl_gyrokinetic_maxwellian_moments *up, 
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range, 
  const struct gkyl_array *fin, struct gkyl_array *density_out)
{
  // compute J*M0 where J is the configurations-space Jacobian
  gkyl_dg_updater_moment_gyrokinetic_advance(up->M0_calc, phase_range, conf_range, 
    fin, up->M0);
  if (up->divide_jacobgeo) {
    // Rescale moment by the inverse of the Jacobian
    gkyl_dg_div_op_range(up->mem, up->conf_basis, 
      0, density_out, 0, up->M0, 0, up->gk_geom->jacobgeo, conf_range);  
  }
  else {
    gkyl_array_set_range(density_out, 1.0, up->M0, conf_range);
  }
}

void 
gkyl_gyrokinetic_maxwellian_moments_advance(struct gkyl_gyrokinetic_maxwellian_moments *up, 
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range, 
  const struct gkyl_array *fin, struct gkyl_array *moms_out)
{
  int vdim = up->vdim;
  int num_conf_basis = up->num_conf_basis;

  // compute J*M0 and J*M1 where J is the configurations-space Jacobian
  gkyl_dg_updater_moment_gyrokinetic_advance(up->M0_calc, phase_range, conf_range, 
    fin, up->M0);
  gkyl_dg_updater_moment_gyrokinetic_advance(up->M1_calc, phase_range, conf_range, 
    fin, up->M1);

  // Isolate u_par by dividing J*M1 by J*M0
  gkyl_dg_div_op_range(up->mem, up->conf_basis, 
    0, up->u_par, 0, up->M1, 0, up->M0, conf_range);

  // Compute J*M1*upar (needed to compute T/m).
  gkyl_dg_mul_op_range(up->conf_basis, 
    0, up->u_par_dot_M1, 0, up->u_par, 0, up->M1, conf_range); 

  // Compute J*M2 = vdim*J*P/m + J*M1*upar.
  gkyl_dg_updater_moment_gyrokinetic_advance(up->M2_calc, phase_range, conf_range, 
    fin, up->pressure);
  // Subtract off J*M1*upar from total J*M2
  gkyl_array_accumulate_range(up->pressure, -1.0, 
    up->u_par_dot_M1, conf_range); 

  // Rescale J*pressure by 1.0/vdim and divide out J*M0 to get T/m, T/m = J*P/(m J*M0). 
  gkyl_array_scale(up->pressure, 1.0/vdim);
  gkyl_dg_div_op_range(up->mem, up->conf_basis, 
    0, up->temperature, 0, up->pressure, 0, up->M0, conf_range);

  if (up->divide_jacobgeo) {
    // Rescale moment by the inverse of the Jacobian and store in moms_out
    gkyl_dg_div_op_range(up->mem, up->conf_basis, 
      0, moms_out, 0, up->M0, 0, up->gk_geom->jacobgeo, conf_range);  
  }
  else {
    gkyl_array_set_range(moms_out, 1.0, up->M0, conf_range);
  }
  // Save the other outputs to moms_out (n, V_drift, T/m):
  gkyl_array_set_offset_range(moms_out, 1.0, up->u_par, 1*num_conf_basis, conf_range);
  gkyl_array_set_offset_range(moms_out, 1.0, up->temperature, 2*num_conf_basis, conf_range);
}

void 
gkyl_gyrokinetic_bi_maxwellian_moments_advance(struct gkyl_gyrokinetic_maxwellian_moments *up, 
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range, 
  const struct gkyl_array *fin, struct gkyl_array *moms_out)
{
  int vdim = up->vdim;
  int num_conf_basis = up->num_conf_basis;

  // compute J*M0 and J*M1 where J is the configurations-space Jacobian
  gkyl_dg_updater_moment_gyrokinetic_advance(up->M0_calc, phase_range, conf_range, 
    fin, up->M0);
  gkyl_dg_updater_moment_gyrokinetic_advance(up->M1_calc, phase_range, conf_range, 
    fin, up->M1);

  // Isolate u_par by dividing J*M1 by J*M0
  gkyl_dg_div_op_range(up->mem, up->conf_basis, 
    0, up->u_par, 0, up->M1, 0, up->M0, conf_range);

  // Compute J*M1*upar (needed to compute T_par/m).
  gkyl_dg_mul_op_range(up->conf_basis, 
    0, up->u_par_dot_M1, 0, up->u_par, 0, up->M1, conf_range); 

  // Compute J*M2_par = J*P_par/m + J*M1*upar.
  gkyl_dg_updater_moment_gyrokinetic_advance(up->M2_par_calc, phase_range, conf_range, 
    fin, up->p_par);
  // Subtract off J*M1*upar from total J*M2_par
  gkyl_array_accumulate_range(up->p_par, -1.0, 
    up->u_par_dot_M1, conf_range); 

  // Compute J*M2_perp = 2*J*P_perp/m.
  gkyl_dg_updater_moment_gyrokinetic_advance(up->M2_perp_calc, phase_range, conf_range, 
    fin, up->p_perp);

  // Rescale J*p_perp by 1.0/2.0 and divide out J*M0 to get T_par/m, T_perp/m from P_par/m, P_perp/m. 
  gkyl_array_scale(up->p_perp, 1.0/2.0);
  gkyl_dg_div_op_range(up->mem, up->conf_basis, 
    0, up->t_par, 0, up->p_par, 0, up->M0, conf_range);
  gkyl_dg_div_op_range(up->mem, up->conf_basis, 
    0, up->t_perp, 0, up->p_perp, 0, up->M0, conf_range);

  if (up->divide_jacobgeo) {
    // Rescale moment by the inverse of the Jacobian and store in moms_out
    gkyl_dg_div_op_range(up->mem, up->conf_basis, 
      0, moms_out, 0, up->M0, 0, up->gk_geom->jacobgeo, conf_range);  
  }
  else {
    gkyl_array_set_range(moms_out, 1.0, up->M0, conf_range);
  }
  // Save the other outputs to moms_out (n, u_par, T_par/m, T_perp/m):
  gkyl_array_set_offset_range(moms_out, 1.0, up->u_par, 1*num_conf_basis, conf_range);
  gkyl_array_set_offset_range(moms_out, 1.0, up->t_par, 2*num_conf_basis, conf_range);
  gkyl_array_set_offset_range(moms_out, 1.0, up->t_perp, 3*num_conf_basis, conf_range);
}

void 
gkyl_gyrokinetic_maxwellian_moments_release(gkyl_gyrokinetic_maxwellian_moments *up)
{
  gkyl_gk_geometry_release(up->gk_geom);
  gkyl_array_release(up->M0);
  gkyl_array_release(up->M1);
  gkyl_array_release(up->u_par);
  gkyl_array_release(up->u_par_dot_M1);
  gkyl_array_release(up->pressure);
  gkyl_array_release(up->temperature);
  gkyl_dg_bin_op_mem_release(up->mem);

  gkyl_dg_updater_moment_gyrokinetic_release(up->M0_calc);
  gkyl_dg_updater_moment_gyrokinetic_release(up->M1_calc);
  gkyl_dg_updater_moment_gyrokinetic_release(up->M2_calc);
  if (up->vdim == 3) {
    gkyl_array_release(up->p_par);
    gkyl_array_release(up->t_par);
    gkyl_array_release(up->p_perp);
    gkyl_array_release(up->t_perp);
    gkyl_dg_updater_moment_gyrokinetic_release(up->M2_par_calc);
    gkyl_dg_updater_moment_gyrokinetic_release(up->M2_perp_calc);
  }

  gkyl_free(up);
}
