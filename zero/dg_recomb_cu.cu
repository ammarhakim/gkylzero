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
#include <gkyl_dg_recomb.h>
#include <gkyl_dg_recomb_priv.h>
#include <gkyl_util.h>
#include <gkyl_const.h>
}

__global__ static void
gkyl_recomb_react_rate_cu_ker(const struct gkyl_dg_recomb *up, const struct gkyl_range conf_rng, const struct gkyl_range adas_rng,
  const struct gkyl_basis *adas_basis, const struct gkyl_dg_prim_vars_type *calc_prim_vars_elc_vtSq,
  const struct gkyl_dg_prim_vars_type *calc_prim_vars_ion_udrift, const struct gkyl_dg_prim_vars_type *calc_prim_vars_ion_vtSq,
  const struct gkyl_array *moms_elc, const struct gkyl_array *moms_ion, struct gkyl_array *vtSq_elc,
  struct gkyl_array *coef_m0, struct gkyl_array *coef_recomb, struct gkyl_array *udrift_ion, struct gkyl_array *vtSq_ion,
  struct gkyl_array *prim_vars_ion, struct gkyl_array *recomb_data, int num_basis,
  enum gkyl_dg_recomb_self type_self, bool all_gk, double mass_elc, double elem_charge, double maxLogTe,
  double minLogTe, double dlogTe, double maxLogM0, double minLogM0, double dlogM0, double resTe, double resM0)
{
  int cidx[GKYL_MAX_CDIM];
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < conf_rng.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&conf_rng, tid, cidx);
    long loc = gkyl_range_idx(&conf_rng, cidx);

    int cdim = conf_rng.ndim;
    const double *moms_elc_d = (const double*) gkyl_array_cfetch(moms_elc, loc);
    const double *m0_elc_d = &moms_elc_d[0];
    
    double *vtSq_elc_d = (double*) gkyl_array_fetch(vtSq_elc, loc);
    double *coef_m0_d = (double*) gkyl_array_fetch(coef_m0, loc);
    double *coef_recomb_d = (double*) gkyl_array_fetch(coef_recomb, loc);

    for (int i=0; i<num_basis; ++i) coef_m0_d[i] = m0_elc_d[i];
    
    calc_prim_vars_elc_vtSq->kernel(calc_prim_vars_elc_vtSq, cidx, moms_elc_d,
				    vtSq_elc_d);

    //Find nearest neighbor for n, Te in ADAS interpolated data
    double cell_av_fac = pow(1.0/sqrt(2.0),cdim);
    double m0_elc_av = m0_elc_d[0]*cell_av_fac;
    double temp_elc_av = vtSq_elc_d[0]*cell_av_fac*mass_elc/elem_charge;
    double log_Te_av = log10(temp_elc_av);
    double log_m0_av = log10(m0_elc_av);
    int m0_idx, t_idx;
    double cell_vals_2d[2];
    double cell_center;
    
    if (log_Te_av < minLogTe) t_idx=1;
    else if (log_Te_av > maxLogTe) t_idx=resTe;
    else t_idx = (log_Te_av - minLogTe)/(dlogTe)+1;
    cell_center = (t_idx - 0.5)*dlogTe + minLogTe;
    cell_vals_2d[0] = 2.0*(log_Te_av - cell_center)/dlogTe; // Te value on cell interval
      
    if (log_m0_av < minLogM0) m0_idx=1;
    else if (log_m0_av > maxLogM0) m0_idx=resM0;
    else m0_idx = (log_m0_av - minLogM0)/(dlogM0)+1;
    cell_center = (m0_idx - 0.5)*dlogM0 + minLogM0;
    cell_vals_2d[1] = 2.0*(log_m0_av - cell_center)/dlogM0; // M0 value on cell interval

    int ad_idx[2] = {t_idx, m0_idx};
    double *recomb_dat_d = (double*) gkyl_array_fetch(recomb_data, gkyl_range_idx(&adas_rng, ad_idx));
    double adas_eval = adas_basis->eval_expand(cell_vals_2d, recomb_dat_d);
    coef_recomb_d[0] = pow(10.0,adas_eval)/cell_av_fac;

    if ((all_gk==false) && (type_self == GKYL_RECOMB_RECVR)) {
      const double *moms_ion_d = (const double*) gkyl_array_cfetch(moms_ion, loc);
      double *udrift_ion_d = (double*) gkyl_array_fetch(udrift_ion, loc);
      double *vtSq_ion_d = (double*) gkyl_array_fetch(vtSq_ion, loc);
      double *prim_vars_ion_d = (double*) gkyl_array_fetch(prim_vars_ion, loc);
      
      // condense the following 2 kernels...
      calc_prim_vars_ion_udrift->kernel(calc_prim_vars_ion_udrift, cidx,
					    moms_ion_d, udrift_ion_d);
      calc_prim_vars_ion_vtSq->kernel(calc_prim_vars_ion_vtSq, cidx,
					    moms_ion_d, vtSq_ion_d);
    }
  }
}

void gkyl_dg_recomb_coll_cu(const struct gkyl_dg_recomb *up,
  const struct gkyl_array *moms_elc, const struct gkyl_array *moms_ion,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, const struct gkyl_array *b_i,
  const struct gkyl_array *f_self, struct gkyl_array *coll_recomb, struct gkyl_array *cflrate)
{

  if ((up->all_gk == false) && (up->type_self == GKYL_RECOMB_RECVR)) {
    // Set auxiliary variable (b_i) for computation of udrift_i
    gkyl_dg_prim_vars_transform_vlasov_gk_set_auxfields(up->calc_prim_vars_ion_udrift, 
      (struct gkyl_dg_prim_vars_auxfields) {.b_i = b_i});
  }
  
  gkyl_recomb_react_rate_cu_ker<<<up->conf_rng->nblocks, up->conf_rng->nthreads>>>(up->on_dev, *up->conf_rng, up->adas_rng, up->basis_on_dev,
    up->calc_prim_vars_elc_vtSq->on_dev, up->calc_prim_vars_ion_udrift->on_dev, up->calc_prim_vars_ion_vtSq->on_dev, 
    moms_elc->on_dev, moms_ion->on_dev, up->vtSq_elc->on_dev, up->coef_m0->on_dev, up->coef_recomb->on_dev,
    up->udrift_ion->on_dev, up->vtSq_ion->on_dev, up->prim_vars_ion->on_dev,
    up->recomb_data->on_dev, up->cbasis->num_basis, up->type_self, up->all_gk, up->mass_elc,
    up->elem_charge, up->maxLogTe, up->minLogTe, up->dlogTe, up->maxLogM0, up->minLogM0,
    up->dlogM0, up->resTe, up->resM0);

  if (up->type_self == GKYL_RECOMB_RECVR) {
    if (up->all_gk) {
      gkyl_proj_gkmaxwellian_on_basis_lab_mom(up->proj_max, up->phase_rng, up->conf_rng, moms_ion, bmag, jacob_tot, up->mass_self, coll_recomb);
    }
    else {
      // Set fmax moments
      gkyl_array_set_offset_range(up->prim_vars_ion, 1., up->udrift_ion, 0, up->conf_rng);
      gkyl_array_set_offset_range(up->prim_vars_ion, 1., up->vtSq_ion, (up->vdim)*up->cbasis->num_basis, up->conf_rng);
      
      // Proj maxwellian on basis
      gkyl_proj_maxwellian_on_basis_prim_mom(up->proj_max, up->phase_rng, up->conf_rng, moms_ion,
					     up->prim_vars_ion, coll_recomb);
    }
  }
  else {
    // copy and weak mult
    gkyl_array_set_range(coll_recomb, -1.0, f_self, up->phase_rng);
  }
  
  gkyl_dg_mul_conf_phase_op_range(up->cbasis, up->pbasis, coll_recomb, up->coef_recomb, coll_recomb,
				  up->conf_rng, up->phase_rng);
  gkyl_dg_mul_conf_phase_op_range(up->cbasis, up->pbasis, coll_recomb, up->coef_m0, coll_recomb,
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