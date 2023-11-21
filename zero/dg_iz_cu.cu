/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_prim_vars_vlasov.h>
#include <gkyl_dg_prim_vars_gyrokinetic.h>
#include <gkyl_dg_prim_vars_transform_vlasov_gk.h>
#include <gkyl_dg_prim_vars_type.h>
#include <gkyl_array_ops.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_iz.h>
#include <gkyl_dg_iz_priv.h>
#include <gkyl_util.h>
#include <gkyl_const.h>
}

__global__ static void
gkyl_iz_react_rate_cu_ker(const struct gkyl_dg_iz *up, const struct gkyl_range conf_rng, const struct gkyl_range adas_rng,
  const struct gkyl_basis adas_basis, const struct gkyl_dg_prim_vars_type *calc_prim_vars_elc_vtSq,
  const struct gkyl_dg_prim_vars_type *calc_prim_vars_donor_gk, const struct gkyl_array *moms_elc,
  const struct gkyl_array *moms_donor, struct gkyl_array *vtSq_elc, struct gkyl_array *prim_vars_donor,
  struct gkyl_array *coef_iz, struct gkyl_array *coef_m0, struct gkyl_array *ioniz_data)
{
  int cidx[GKYL_MAX_CDIM];
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < conf_rng.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&conf_rng, tid, cidx);
    long loc = gkyl_range_idx(&conf_rng, cidx);

    const double *moms_elc_d = (const double*) gkyl_array_cfetch(moms_elc, loc);
    const double *m0_elc_d = &moms_elc_d[0];
    
    double *vtSq_elc_d = (double*) gkyl_array_fetch(vtSq_elc, loc);
    double *coef_m0_d = (double*) gkyl_array_fetch(coef_m0, loc);
    double *coef_iz_d = (double*) gkyl_array_fetch(coef_iz, loc);

    for (int i=0; i<up->cbasis->num_basis; ++i) coef_m0_d[i] = m0_elc_d[i];
    
    calc_prim_vars_elc_vtSq->kernel(calc_prim_vars_elc_vtSq, cidx, moms_elc_d,
				    vtSq_elc_d);

    if ((up->type_self == GKYL_IZ_ELC) || (up->type_self == GKYL_IZ_ION)) {
      const double *moms_donor_d = (const double*) gkyl_array_cfetch(moms_donor, loc);
      const double *m0_donor_d = &moms_donor_d[0];
      double *prim_vars_donor_d = (double*) gkyl_array_fetch(prim_vars_donor, loc);
      calc_prim_vars_donor_gk->kernel(calc_prim_vars_donor_gk, cidx,
				      moms_donor_d, prim_vars_donor_d);
      if (up->type_self == GKYL_IZ_ELC) {
      	for (int i=0; i<up->cbasis->num_basis; ++i) coef_m0_d[i] = m0_donor_d[i];
      }
    }   
    
    //Find nearest neighbor for n, Te in ADAS interpolated data
    double cell_av_fac = pow(1.0/sqrt(2.0),up->cdim);
    double m0_elc_av = m0_elc_d[0]*cell_av_fac;
    double temp_elc_av = vtSq_elc_d[0]*cell_av_fac*up->mass_elc/up->elem_charge;
    double log_Te_av = log10(temp_elc_av);
    double log_m0_av = log10(m0_elc_av);
    double cell_val_t;
    double cell_val_m0;
    int m0_idx, t_idx;
    double cell_vals_2d[2];
    double cell_center;
    
    if (log_Te_av < up->minLogTe) t_idx=1;
    else if (log_Te_av > up->maxLogTe) t_idx=up->resTe;
    else t_idx = (log_Te_av - up->minLogTe)/(up->dlogTe)+1;
    cell_center = (t_idx - 0.5)*up->dlogTe + up->minLogTe;
    cell_vals_2d[0] = 2.0*(log_Te_av - cell_center)/up->dlogTe; // Te value on cell interval
      
    if (log_m0_av < up->minLogM0) m0_idx=1;
    else if (log_m0_av > up->maxLogM0) m0_idx=up->resM0;
    else m0_idx = (log_m0_av - up->minLogM0)/(up->dlogM0)+1;
    cell_center = (m0_idx - 0.5)*up->dlogM0 + up->minLogM0;
    cell_vals_2d[1] = 2.0*(log_m0_av - cell_center)/up->dlogM0; // M0 value on cell interval
 
    if ((up->E/temp_elc_av >= 3./2.) || (m0_elc_av <= 0.) || (temp_elc_av <= 0.)) {
      coef_iz_d[0] = 0.0;
    }
    else {
      double *iz_dat_d = (double*) gkyl_array_fetch(up->ioniz_data, gkyl_range_idx(&up->adas_rng, (int[2]) {t_idx,m0_idx}));
      double adas_eval = up->adas_basis.eval_expand(cell_vals_2d, iz_dat_d);
      coef_iz_d[0] = pow(10.0,adas_eval)/cell_av_fac;
    }
  }
}

void gkyl_dg_iz_coll_cu(const struct gkyl_dg_iz *up,
  const struct gkyl_array *moms_elc, const struct gkyl_array *moms_donor,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, const struct gkyl_array *b_i,
  const struct gkyl_array *f_self, struct gkyl_array *coll_iz, struct gkyl_array *cflrate)
{
  if ((up->all_gk==false) && ((up->type_self == GKYL_IZ_ELC) || (up->type_self == GKYL_IZ_ION))) {
    // Set auxiliary variable (b_i) for computation of gk neut prim vars
    gkyl_dg_prim_vars_transform_vlasov_gk_set_auxfields(up->calc_prim_vars_donor, 
      (struct gkyl_dg_prim_vars_auxfields) {.b_i = b_i});
  }x

  gkyl_iz_react_rate_cu_ker<<<up->conf_rng->nblocks, up->conf_rng->nthreads>>>(up->on_dev, *up->conf_rng, up->adas_rng,
    up->adas_basis, up->calc_prim_vars_elc_vtSq->on_dev, up->calc_prim_vars_donor->on_dev, 
    moms_elc->on_dev, moms_donor->on_dev, up->vtSq_elc->on_dev, up->prim_vars_donor->on_dev,
    up->coef_iz->on_dev, up->coef_m0->on_dev, up->ioniz_data->on_dev);

  if (up->type_self == GKYL_IZ_ELC) {
     
    // Calculate vt_sq_iz
    gkyl_array_copy_range(up->vtSq_iz, up->vtSq_elc, up->conf_rng);
    gkyl_array_scale_range(up->vtSq_iz, 1/2.0, up->conf_rng);
    gkyl_array_shiftc(up->vtSq_iz, -up->E*up->elem_charge/(3*up->mass_elc)*pow(sqrt(2),up->cdim), 0);

    // Set fmax moments
    gkyl_array_set_offset_range(up->prim_vars_fmax, 1., up->prim_vars_donor, 0, up->conf_rng);
    gkyl_array_set_offset_range(up->prim_vars_fmax, 1., up->vtSq_iz, up->cbasis->num_basis, up->conf_rng);

    // Proj maxwellian on basis
    gkyl_proj_gkmaxwellian_on_basis_prim_mom(up->proj_max, up->phase_rng, up->conf_rng, moms_elc,
    					     up->prim_vars_fmax, bmag, jacob_tot, up->mass_elc, coll_iz);
    
    // copy, scale and accumulate
    gkyl_array_scale_range(coll_iz, 2.0, up->phase_rng);
    gkyl_array_accumulate_range(coll_iz, -1.0, f_self, up->phase_rng);

  }
  else if (up->type_self == GKYL_IZ_ION) {
    // Proj maxwellian on basis (doesn't assume same phase grid, even if GK)
    gkyl_proj_gkmaxwellian_on_basis_prim_mom(up->proj_max, up->phase_rng, up->conf_rng, moms_donor,
    					     up->prim_vars_donor, bmag, jacob_tot, up->mass_ion, coll_iz);
  }
  else if (up->type_self == GKYL_IZ_DONOR) {
    // neut coll_iz = -f_n
    gkyl_array_set_range(coll_iz, -1.0, f_self, up->phase_rng);
  }

  // coll_iz = coef_iz*coll_iz
  gkyl_dg_mul_conf_phase_op_range(up->cbasis, up->pbasis, coll_iz, up->coef_iz, coll_iz,
  				    up->conf_rng, up->phase_rng);
  gkyl_dg_mul_conf_phase_op_range(up->cbasis, up->pbasis, coll_iz, up->coef_m0, coll_iz,
  				    up->conf_rng, up->phase_rng);
  
  
  // cfl calculation
  //struct gkyl_range vel_rng;
  /* gkyl_range_deflate(&vel_rng, up->phase_rng, rem_dir, conf_iter.idx); */
  /* gkyl_range_iter_no_split_init(&vel_iter, &vel_rng); */
  /* // cfl associated with reaction is a *phase space* cfl */
  /* // Need to loop over velocity space for each configuration space cell */
  /* // to get total cfl rate in each phase space cell */
  /* while (gkyl_range_iter_next(&vel_iter)) { */
  /*   long cfl_idx = gkyl_range_idx(&vel_rng, vel_iter.idx); */
  /*   double *cflrate_d = gkyl_array_fetch(cflrate, cfl_idx); */
  /*   cflrate_d[0] += cflr; // frequencies are additive */
  /* } */
}
