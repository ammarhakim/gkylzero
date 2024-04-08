#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_dg_updater_moment.h>
#include <gkyl_vlasov_lte_moments.h>
#include <gkyl_vlasov_lte_moments_priv.h>
#include <gkyl_mom_vlasov_sr.h>

struct gkyl_vlasov_lte_moments*
gkyl_vlasov_lte_moments_inew(const struct gkyl_vlasov_lte_moments_inp *inp)
{
  gkyl_vlasov_lte_moments *up = gkyl_malloc(sizeof(*up));

  up->conf_basis = *inp->conf_basis;
  up->phase_basis = *inp->phase_basis;
  up->num_conf_basis = inp->conf_basis->num_basis;
  up->vdim = up->phase_basis.ndim - up->conf_basis.ndim;
  up->model_id = inp->model_id;
  up->mass = inp->mass;

  long conf_local_ncells = inp->conf_range->volume;
  long conf_local_ext_ncells = inp->conf_range_ext->volume;

  if (inp->use_gpu) {
    up->M0 = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
    up->M1i = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->vdim*up->num_conf_basis, conf_local_ext_ncells);
    up->V_drift = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->vdim*up->num_conf_basis, conf_local_ext_ncells);
    up->V_drift_dot_M1i = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
    up->pressure = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
    up->temperature = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
    up->mem = gkyl_dg_bin_op_mem_cu_dev_new(conf_local_ncells, up->num_conf_basis);
  }
  else {
    up->M0 = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
    up->M1i = gkyl_array_new(GKYL_DOUBLE, up->vdim*up->num_conf_basis, conf_local_ext_ncells);
    up->V_drift = gkyl_array_new(GKYL_DOUBLE, up->vdim*up->num_conf_basis, conf_local_ext_ncells);
    up->V_drift_dot_M1i = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
    up->pressure = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
    up->temperature = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
    up->mem = gkyl_dg_bin_op_mem_new(conf_local_ncells, up->num_conf_basis);
  }

  if (up->model_id == GKYL_MODEL_SR) {
    if (inp->use_gpu) {
      up->GammaV2 = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
      up->GammaV_inv = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
      up->M0_minus_V_drift_dot_M1i = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
    }
    else {
      up->GammaV2 = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
      up->GammaV_inv = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
      up->M0_minus_V_drift_dot_M1i = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
    }
    // Set auxiliary fields for moment updates. 
    struct gkyl_mom_vlasov_sr_auxfields sr_inp = {.p_over_gamma = inp->p_over_gamma, 
      .gamma = inp->gamma, .gamma_inv = inp->gamma_inv, .V_drift = up->V_drift, 
      .GammaV2 = up->GammaV2, .GammaV_inv = up->GammaV_inv};  
    // Moment calculator for needed moments (M0, M1i, and P for relativistic)
    up->M0_calc = gkyl_dg_updater_moment_new(inp->phase_grid, inp->conf_basis,
      inp->phase_basis, inp->conf_range, inp->vel_range, up->model_id, &sr_inp, "M0", 0, up->mass, inp->use_gpu);
    up->M1i_calc = gkyl_dg_updater_moment_new(inp->phase_grid, inp->conf_basis,
      inp->phase_basis, inp->conf_range, inp->vel_range, up->model_id, &sr_inp, "M1i", 0, up->mass, inp->use_gpu);
    up->Pcalc = gkyl_dg_updater_moment_new(inp->phase_grid, inp->conf_basis,
      inp->phase_basis, inp->conf_range, inp->vel_range, up->model_id, &sr_inp, "Pressure", 0, up->mass, inp->use_gpu);
  }
  else {
    // Moment calculator for needed moments (M0, M1i, and M2 for non-relativistic)
    // Note: auxiliary field input is NULL (not used by non-relativistic simulations)
    up->M0_calc = gkyl_dg_updater_moment_new(inp->phase_grid, inp->conf_basis,
      inp->phase_basis, inp->conf_range, inp->vel_range, up->model_id, 0, "M0", 0, up->mass, inp->use_gpu);
    up->M1i_calc = gkyl_dg_updater_moment_new(inp->phase_grid, inp->conf_basis,
      inp->phase_basis, inp->conf_range, inp->vel_range, up->model_id, 0, "M1i", 0, up->mass, inp->use_gpu);
    up->Pcalc = gkyl_dg_updater_moment_new(inp->phase_grid, inp->conf_basis,
      inp->phase_basis, inp->conf_range, inp->vel_range, up->model_id, 0, "M2", 0, up->mass, inp->use_gpu);    
  }
  return up;
}

void 
gkyl_vlasov_lte_density_moment_advance(struct gkyl_vlasov_lte_moments *lte_moms, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local, 
  const struct gkyl_array *fin, struct gkyl_array *density_out)
{
  int vdim = lte_moms->vdim;
  // compute lab frame moment M0
  gkyl_dg_updater_moment_advance(lte_moms->M0_calc, phase_local, conf_local, 
    fin, lte_moms->M0);

  // If we are relativistic, we need to compute the relevant Lorentz 
  // boost factors and perform the Lorentz transformation to go from
  // the lab frame to the stationary frame. 
  if (lte_moms->model_id == GKYL_MODEL_SR) {
    // Need V_drift in relativity to compute the Lorentz boost factors
    gkyl_dg_updater_moment_advance(lte_moms->M1i_calc, phase_local, conf_local, 
      fin, lte_moms->M1i);
    // Isolate drift velocity by dividing M1i by M0
    for (int d = 0; d < vdim; ++d) {
      gkyl_dg_div_op_range(lte_moms->mem, lte_moms->conf_basis, 
        d, lte_moms->V_drift,
        d, lte_moms->M1i, 0, lte_moms->M0, conf_local);
    }
    // Compute V_drift dot M1i (needed to compute stationary frame moments).
    gkyl_array_clear(lte_moms->V_drift_dot_M1i, 0.0);
    gkyl_dg_dot_product_op_range(lte_moms->conf_basis, 
      lte_moms->V_drift_dot_M1i, lte_moms->V_drift, lte_moms->M1i, conf_local); 
        
    gkyl_calc_sr_vars_GammaV_inv(&lte_moms->conf_basis, &lte_moms->phase_basis,
      conf_local, lte_moms->V_drift, lte_moms->GammaV_inv);

    // ( n = GammaV*(M0 - V_drift dot M1i) ) Lorentz transform to our fluid-stationary density 
    // This expression follows from the fact that M0 = GammaV*n and M1i = GammaV*n*V_drift so
    // n = GammaV^2*n*(1 - V_drift^2) = n*(1 - V_drift^2)/(1 - V_drift^2) = n
    gkyl_array_set(lte_moms->M0_minus_V_drift_dot_M1i, 1.0, lte_moms->M0);
    gkyl_array_accumulate_range(lte_moms->M0_minus_V_drift_dot_M1i, -1.0, 
      lte_moms->V_drift_dot_M1i, conf_local);
    gkyl_dg_div_op_range(lte_moms->mem,lte_moms->conf_basis, 0, density_out, 
      0, lte_moms->M0_minus_V_drift_dot_M1i, 0, lte_moms->GammaV_inv, conf_local);
  }
  else {
    gkyl_array_set_range(density_out, 1.0, lte_moms->M0, conf_local);
  }
}

void 
gkyl_vlasov_lte_moments_advance(struct gkyl_vlasov_lte_moments *lte_moms, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local, 
  const struct gkyl_array *fin, struct gkyl_array *moms_out)
{
  int vdim = lte_moms->vdim;
  int num_conf_basis = lte_moms->num_conf_basis;

  // compute lab frame moments M0 and M1i
  gkyl_dg_updater_moment_advance(lte_moms->M0_calc, phase_local, conf_local, 
    fin, lte_moms->M0);
  gkyl_dg_updater_moment_advance(lte_moms->M1i_calc, phase_local, conf_local, 
    fin, lte_moms->M1i);

  // Isolate drift velocity by dividing M1i by M0
  for (int d = 0; d < vdim; ++d) {
    gkyl_dg_div_op_range(lte_moms->mem, lte_moms->conf_basis, 
      d, lte_moms->V_drift,
      d, lte_moms->M1i, 0, lte_moms->M0, conf_local);
  }
  // Compute V_drift dot M1i (needed to compute stationary frame moments).
  gkyl_array_clear(lte_moms->V_drift_dot_M1i, 0.0);
  gkyl_dg_dot_product_op_range(lte_moms->conf_basis, 
    lte_moms->V_drift_dot_M1i, lte_moms->V_drift, lte_moms->M1i, conf_local); 

  // If we are relativistic, we need to compute the relevant Lorentz 
  // boost factors and perform the Lorentz transformation to go from
  // the lab frame to the stationary frame. 
  if (lte_moms->model_id == GKYL_MODEL_SR) {
    // (GammaV^2 = 1/(1-V_drift^2)) 
    gkyl_calc_sr_vars_GammaV2(&lte_moms->conf_basis, &lte_moms->phase_basis,
      conf_local, lte_moms->V_drift, lte_moms->GammaV2);
    // (GammaV_inv = sqrt(1-V_drift^2))
    gkyl_calc_sr_vars_GammaV_inv(&lte_moms->conf_basis, &lte_moms->phase_basis,
      conf_local, lte_moms->V_drift, lte_moms->GammaV_inv);

    // Compute the pressure moment.
    // This moment is computed *in the stationary frame* in the relativistic moment calculator.
    // We find computing this moment *in the stationary frame* to be more accurate than 
    // computing the lab frame moment and then Lorentz transforming to the stationary frame. 
    gkyl_dg_updater_moment_advance(lte_moms->Pcalc, phase_local, conf_local, 
      fin, lte_moms->pressure);

    // ( n = GammaV*(M0 - V_drift dot M1i) ) Lorentz transform to our fluid-stationary density 
    // This expression follows from the fact that M0 = GammaV*n and M1i = GammaV*n*V_drift so
    // n = GammaV^2*n*(1 - V_drift^2) = n*(1 - V_drift^2)/(1 - V_drift^2) = n
    gkyl_array_set(lte_moms->M0_minus_V_drift_dot_M1i, 1.0, lte_moms->M0);
    gkyl_array_accumulate_range(lte_moms->M0_minus_V_drift_dot_M1i, -1.0, 
      lte_moms->V_drift_dot_M1i, conf_local);
    gkyl_dg_div_op_range(lte_moms->mem,lte_moms->conf_basis, 0, moms_out, 
      0, lte_moms->M0_minus_V_drift_dot_M1i, 0, lte_moms->GammaV_inv, conf_local);
  }
  else {
    // Compute the lab frame M2 = vdim*P/m + V_drift dot M1i.
    gkyl_dg_updater_moment_advance(lte_moms->Pcalc, phase_local, conf_local, 
      fin, lte_moms->pressure);
    // Subtract off V_drift dot M1i from total M2
    gkyl_array_accumulate_range(lte_moms->pressure, -1.0, 
      lte_moms->V_drift_dot_M1i, conf_local); 

    // Rescale pressure by 1.0/vdim and set the first component of moms_out to be the density. 
    gkyl_array_scale(lte_moms->pressure, 1.0/vdim);
    gkyl_array_set_range(moms_out, 1.0, lte_moms->M0, conf_local);
  }
  // ( T/m = P/(mn) ) 
  gkyl_dg_div_op_range(lte_moms->mem, lte_moms->conf_basis, 
    0, lte_moms->temperature,
    0, lte_moms->pressure, 0, moms_out, conf_local);
  // Save the outputs to moms_out (n, V_drift, T/m):
  gkyl_array_set_offset_range(moms_out, 1.0, lte_moms->V_drift, 1*num_conf_basis, conf_local);
  gkyl_array_set_offset_range(moms_out, 1.0, lte_moms->temperature, (vdim+1)*num_conf_basis, conf_local);
}

void 
gkyl_vlasov_lte_moments_release(gkyl_vlasov_lte_moments *lte_moms)
{
  gkyl_array_release(lte_moms->M0);
  gkyl_array_release(lte_moms->M1i);
  gkyl_array_release(lte_moms->V_drift);
  gkyl_array_release(lte_moms->V_drift_dot_M1i);
  gkyl_array_release(lte_moms->pressure);
  gkyl_array_release(lte_moms->temperature);
  gkyl_dg_bin_op_mem_release(lte_moms->mem);
  if (lte_moms->model_id == GKYL_MODEL_SR) {
    gkyl_array_release(lte_moms->M0_minus_V_drift_dot_M1i);
    gkyl_array_release(lte_moms->GammaV_inv);
    gkyl_array_release(lte_moms->GammaV2);    
  }
  gkyl_dg_updater_moment_release(lte_moms->M0_calc);
  gkyl_dg_updater_moment_release(lte_moms->M1i_calc);
  gkyl_dg_updater_moment_release(lte_moms->Pcalc);

  gkyl_free(lte_moms);
}
