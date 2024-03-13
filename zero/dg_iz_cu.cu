/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_prim_vars_vlasov.h>
#include <gkyl_dg_prim_vars_gyrokinetic.h>
#include <gkyl_dg_prim_vars_transform.h>
#include <gkyl_dg_prim_vars_type.h>
#include <gkyl_array_ops.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_iz.h>
#include <gkyl_dg_iz_priv.h>
#include <gkyl_util.h>
#include <gkyl_const.h>
#include <gkyl_array_ops_priv.h>
}

__global__ static void
gkyl_iz_react_rate_cu_ker(const struct gkyl_dg_iz *up, const struct gkyl_range conf_rng, const struct gkyl_range adas_rng,
  const struct gkyl_basis *adas_basis, const struct gkyl_dg_prim_vars_type *calc_prim_vars_elc_vtSq,
  const struct gkyl_dg_prim_vars_type *calc_prim_vars_donor_gk, const struct gkyl_array* moms_elc,
  const struct gkyl_array* moms_donor, struct gkyl_array* vtSq_elc, struct gkyl_array* vtSq_iz, struct gkyl_array* prim_vars_donor,
  struct gkyl_array* coef_iz, struct gkyl_array* ioniz_data, int num_basis, enum gkyl_react_self_type type_self,
  double mass_elc, double elem_charge, double E, double maxLogTe, double minLogTe, double dlogTe,
  double maxLogM0, double minLogM0, double dlogM0, int resTe, int resM0, long nc)
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
    double *vtSq_iz_d = (double*) gkyl_array_fetch(vtSq_iz, loc);
    double *coef_iz_d = (double*) gkyl_array_fetch(coef_iz, loc);
    
    calc_prim_vars_elc_vtSq->kernel(calc_prim_vars_elc_vtSq, cidx, moms_elc_d,
				    vtSq_elc_d);

    if ((type_self == GKYL_SELF_ELC) || (type_self == GKYL_SELF_ION)) {
      const double *moms_donor_d = (const double*) gkyl_array_cfetch(moms_donor, loc);
      const double *m0_donor_d = &moms_donor_d[0];
      double *prim_vars_donor_d = (double*) gkyl_array_fetch(prim_vars_donor, loc);
      calc_prim_vars_donor_gk->kernel(calc_prim_vars_donor_gk, cidx,
				      moms_donor_d, prim_vars_donor_d);
    }   
    
    //Find nearest neighbor for n, Te in ADAS interpolated data
    double cell_av_fac = pow(1.0/sqrt(2.0),cdim);
    double m0_elc_av = m0_elc_d[0]*cell_av_fac;
    double temp_elc_av = vtSq_elc_d[0]*cell_av_fac*mass_elc/elem_charge;
    double log_Te_av = log10(temp_elc_av);
    double log_m0_av = log10(m0_elc_av);
    double cell_val_t;
    double cell_val_m0;
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
 
    if ((E/temp_elc_av >= 3./2.) || (m0_elc_av <= 0.)) {
      coef_iz_d[0] = 0.0;
    }
    else {
      int ad_idx[2] = {t_idx, m0_idx};
      double *iz_dat_d = (double*) gkyl_array_fetch(ioniz_data, gkyl_range_idx(&adas_rng, ad_idx));
      double adas_eval = adas_basis->eval_expand(cell_vals_2d, iz_dat_d);
      coef_iz_d[0] = pow(10.0,adas_eval)/cell_av_fac;
      if (temp_elc_av <= 0.) {
	array_set1(nc, vtSq_iz_d, 0.5, vtSq_elc_d);
      	vtSq_iz_d[0] = vtSq_iz_d[0] - E*elem_charge/(3*mass_elc*cell_av_fac);
      }
      else {
	vtSq_iz_d[0] = E/(6.0*mass_elc)/cell_av_fac;
      }
    }
  }
}

void gkyl_dg_iz_coll_cu(const struct gkyl_dg_iz *up, const struct gkyl_array *moms_elc,
  const struct gkyl_array *moms_donor, const struct gkyl_array *b_i,
  struct gkyl_array *vtSq_iz, struct gkyl_array *prim_vars_donor,		 
  struct gkyl_array *coef_iz, struct gkyl_array *cflrate)
{
  if ((up->all_gk==false) && ((up->type_self == GKYL_SELF_ELC) || (up->type_self == GKYL_SELF_ION))) {
    // Set auxiliary variable (b_i) for computation of gk neut prim vars
    gkyl_dg_prim_vars_transform_set_auxfields(up->calc_prim_vars_donor, 
      (struct gkyl_dg_prim_vars_auxfields) {.b_i = b_i});
  }

  long nc = vtSq_iz->ncomp; 
  gkyl_iz_react_rate_cu_ker<<<up->conf_rng->nblocks, up->conf_rng->nthreads>>>(up->on_dev, *up->conf_rng, up->adas_rng,
    up->basis_on_dev, up->calc_prim_vars_elc_vtSq->on_dev, up->calc_prim_vars_donor->on_dev, 
    moms_elc->on_dev, moms_donor->on_dev, up->vtSq_elc->on_dev, vtSq_iz->on_dev, prim_vars_donor->on_dev,
    coef_iz->on_dev, up->ioniz_data->on_dev, up->cbasis->num_basis,
    up->type_self, up->mass_elc, up->elem_charge, up->E, up->maxLogTe, up->minLogTe,
    up->dlogTe, up->maxLogM0, up->minLogM0, up->dlogM0, up->resTe, up->resM0, nc);
  
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
