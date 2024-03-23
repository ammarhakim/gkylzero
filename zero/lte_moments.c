#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_dg_updater_moment.h>
#include <gkyl_lte_moments.h>
#include <gkyl_lte_moments_priv.h>
#include <gkyl_mom_vlasov_sr.h>

struct gkyl_lte_moments*
gkyl_lte_moments_inew(const struct gkyl_lte_moments_vlasov_inp *inp)
{
  gkyl_lte_moments *up = gkyl_malloc(sizeof(*up));

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
      up->Gamma = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
      up->GammaV2 = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
      up->Gamma_inv = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
      up->M0_minus_V_drift_dot_M1i = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
    }
    else {
      up->Gamma = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
      up->GammaV2 = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
      up->Gamma_inv = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
      up->M0_minus_V_drift_dot_M1i = gkyl_array_new(GKYL_DOUBLE, up->num_conf_basis, conf_local_ext_ncells);
    }
    // Set auxiliary fields for moment updates. 
    struct gkyl_mom_vlasov_sr_auxfields sr_inp = {.p_over_gamma = inp->p_over_gamma, 
      .gamma = inp->gamma, .gamma_inv = inp->gamma_inv, .V_drift = up->V_drift, 
      .GammaV2 = up->GammaV2, .GammaV_inv = up->Gamma_inv};  
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
gkyl_lte_moments_advance(struct gkyl_lte_moments *maxwell_moms, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local, 
  const struct gkyl_array *fin, struct gkyl_array *moms)
{
  int vdim = maxwell_moms->vdim;
  int num_conf_basis = maxwell_moms->num_conf_basis;

  // compute lab frame moments M0 and M1i
  gkyl_dg_updater_moment_advance(maxwell_moms->M0_calc, phase_local, conf_local, 
    fin, maxwell_moms->M0);
  gkyl_dg_updater_moment_advance(maxwell_moms->M1i_calc, phase_local, conf_local, 
    fin, maxwell_moms->M1i);

  // Isolate drift velocity by dividing M1i by M0
  for (int d = 0; d < vdim; ++d) {
    gkyl_dg_div_op_range(maxwell_moms->mem, maxwell_moms->conf_basis, 
      d, maxwell_moms->V_drift,
      d, maxwell_moms->M1i, 0, maxwell_moms->M0, conf_local);
  }
  // Compute V_drift dot M1i (needed to compute stationary frame moments).
  gkyl_array_clear(maxwell_moms->V_drift_dot_M1i, 0.0);
  gkyl_dg_dot_product_op_range(maxwell_moms->conf_basis, 
    maxwell_moms->V_drift_dot_M1i, maxwell_moms->V_drift, maxwell_moms->M1i, conf_local); 

  // If we are relativistic, we need to compute the relevant Lorentz 
  // boost factors and perform the Lorentz transformation to go from
  // the lab frame to the stationary frame. 
  if (maxwell_moms->model_id == GKYL_MODEL_SR) {
    // (Gamma = 1/sqrt(1-V_drift^2)) 
    gkyl_calc_sr_vars_Gamma(&maxwell_moms->conf_basis, &maxwell_moms->phase_basis,
         conf_local, maxwell_moms->V_drift, maxwell_moms->Gamma);
    // (Gamma^2 = 1/(1-V_drift^2)) 
    gkyl_calc_sr_vars_Gamma2(&maxwell_moms->conf_basis, &maxwell_moms->phase_basis,
      conf_local, maxwell_moms->V_drift, maxwell_moms->GammaV2);
    // (Gamma_inv = sqrt(1-V_drift^2))
    gkyl_calc_sr_vars_Gamma_inv(&maxwell_moms->conf_basis, &maxwell_moms->phase_basis,
      conf_local, maxwell_moms->V_drift, maxwell_moms->Gamma_inv);

    // Compute the pressure moment.
    // This moment is computed *in the stationary frame* in the relativistic moment calculator.
    // We find computing this moment *in the stationary frame* to be more accurate than 
    // computing the lab frame moment and then Lorentz transforming to the stationary frame. 
    gkyl_dg_updater_moment_advance(maxwell_moms->Pcalc, phase_local, conf_local, 
      fin, maxwell_moms->pressure);

    // ( n = Gamma*(M0 - V_drift dot M1i) ) Lorentz transform to our fluid-stationary density 
    // This expression follows from the fact that M0 = Gamma*n and M1i = Gamma*n*V_drift so
    // n = Gamma^2*n*(1 - V_drift^2) = n*(1 - V_drift^2)/(1 - V_drift^2) = n
    gkyl_array_set(maxwell_moms->M0_minus_V_drift_dot_M1i, 1.0, maxwell_moms->M0);
    gkyl_array_accumulate_range(maxwell_moms->M0_minus_V_drift_dot_M1i, -1.0, 
      maxwell_moms->V_drift_dot_M1i, conf_local);
    //gkyl_dg_mul_op_range(maxwell_moms->conf_basis, 0, moms, 
    //  0, maxwell_moms->Gamma, 0, maxwell_moms->M0_minus_V_drift_dot_M1i, conf_local);
    gkyl_dg_div_op_range(maxwell_moms->mem,maxwell_moms->conf_basis, 0, moms, 
      0, maxwell_moms->M0_minus_V_drift_dot_M1i, 0, maxwell_moms->Gamma_inv, conf_local);
  }
  else {
    // Compute the lab frame M2 = vdim*P/m + V_drift dot M1i.
    gkyl_dg_updater_moment_advance(maxwell_moms->Pcalc, phase_local, conf_local, 
      fin, maxwell_moms->pressure);
    // Subtract off V_drift dot M1i from total M2
    gkyl_array_accumulate_range(maxwell_moms->pressure, -1.0, 
      maxwell_moms->V_drift_dot_M1i, conf_local); 

    // Rescale pressure by 1.0/vdim and set the first component of moms to be the density. 
    gkyl_array_scale(maxwell_moms->pressure, 1.0/vdim);
    gkyl_array_set_range(moms, 1.0, maxwell_moms->M0, conf_local);
  }
  // ( T/m = P/(mn) ) 
  gkyl_dg_div_op_range(maxwell_moms->mem, maxwell_moms->conf_basis, 
    0, maxwell_moms->temperature,
    0, maxwell_moms->pressure, 0, moms, conf_local);
  // Save the outputs to moms (n, V_drift, T/m):
  gkyl_array_set_offset_range(moms, 1.0, maxwell_moms->V_drift, 1*num_conf_basis, conf_local);
  gkyl_array_set_offset_range(moms, 1.0, maxwell_moms->temperature, (vdim+1)*num_conf_basis, conf_local);
}

void 
gkyl_lte_moments_release(gkyl_lte_moments *maxwell_moms)
{
  gkyl_array_release(maxwell_moms->M0);
  gkyl_array_release(maxwell_moms->M1i);
  gkyl_array_release(maxwell_moms->V_drift);
  gkyl_array_release(maxwell_moms->V_drift_dot_M1i);
  gkyl_array_release(maxwell_moms->pressure);
  gkyl_array_release(maxwell_moms->temperature);
  gkyl_dg_bin_op_mem_release(maxwell_moms->mem);
  if (maxwell_moms->model_id == GKYL_MODEL_SR) {
    gkyl_array_release(maxwell_moms->M0_minus_V_drift_dot_M1i);
    gkyl_array_release(maxwell_moms->Gamma_inv);
    gkyl_array_release(maxwell_moms->Gamma);
    gkyl_array_release(maxwell_moms->GammaV2);    
  }
  gkyl_dg_updater_moment_release(maxwell_moms->M0_calc);
  gkyl_dg_updater_moment_release(maxwell_moms->M1i_calc);
  gkyl_dg_updater_moment_release(maxwell_moms->Pcalc);

  gkyl_free(maxwell_moms);
}
