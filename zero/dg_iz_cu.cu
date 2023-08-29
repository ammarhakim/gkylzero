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
}

__global__ static void
gkyl_iz_react_rate_cu_ker(const struct gkyl_dg_iz *up, const struct gkyl_range conf_rng,
  const struct gkyl_dg_prim_vars_type *calc_prim_vars_elc_vtSq, const struct gkyl_dg_prim_vars_type *calc_prim_vars_neut_gk,
  const struct gkyl_array *moms_elc, const struct gkyl_array *moms_neut,
  struct gkyl_array *vtSq_elc, struct gkyl_array *prim_vars_neut, struct gkyl_array *coef_iz, 
  struct gkyl_array *M0q, struct gkyl_array *Teq, struct gkyl_array *ioniz_data)
{
  int cidx[GKYL_MAX_CDIM];
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < conf_rng.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&conf_rng, tid, cidx);
    long loc = gkyl_range_idx(&conf_rng, cidx);

    const double *moms_elc_d = (const double*) gkyl_array_cfetch(moms_elc, loc);
    const double *moms_neut_d = (const double*) gkyl_array_cfetch(moms_neut, loc);
    const double *m0_elc_d = &moms_elc_d[0];

    double *vtSq_elc_d = (double*) gkyl_array_fetch(vtSq_elc, loc);
    double *prim_vars_neut_d = (double*) gkyl_array_fetch(prim_vars_neut, loc);
    double *coef_iz_d = (double*) gkyl_array_fetch(coef_iz, loc); 
    calc_prim_vars_elc_vtSq->kernel(calc_prim_vars_elc_vtSq, cidx, moms_elc_d, vtSq_elc_d);
    calc_prim_vars_neut_gk->kernel(calc_prim_vars_neut_gk, cidx, moms_neut_d, prim_vars_neut_d);
    //Find nearest neighbor for n, Te in ADAS interpolated data
    double cell_av_fac = pow(1.0/sqrt(2.0),up->cdim);
    double m0_elc_av = m0_elc_d[0]*cell_av_fac;
    double temp_elc_av = vtSq_elc_d[0]*cell_av_fac*up->mass_elc/up->elem_charge;
    double diff1 = 0;
    double diff2 = 0;
    int m0_idx, t_idx, q_idx;
    int qtotal = up->resM0*up->resTe;
    
    if (log10(m0_elc_av) < up->minLogM0) m0_idx=0;
    else if (log10(m0_elc_av) > up->maxLogM0) m0_idx=qtotal-1;
    else {
      m0_idx = (log10(m0_elc_av) - up->minLogM0)/(up->dlogM0);
      double *M0q_1 = (double*) gkyl_array_fetch(M0q, m0_idx*up->resTe);
      double *M0q_2 = (double*) gkyl_array_fetch(M0q, (1+m0_idx)*up->resTe);
      diff1 = fabs(M0q_1[0]-log10(m0_elc_av));
      diff2 = fabs(M0q_2[0]-log10(m0_elc_av));
      if (diff1 > diff2) m0_idx += 1;
    }

    if (log10(temp_elc_av) < up->minLogTe) t_idx=0;
    else if (log10(temp_elc_av) > up->maxLogTe) t_idx=up->resTe-1;
    else {
      t_idx = (log10(temp_elc_av) - up->minLogTe)/(up->dlogTe);
      double *Teq_1 = (double*) gkyl_array_fetch(Teq, t_idx);
      double *Teq_2 = (double*) gkyl_array_fetch(Teq, t_idx+1); 
      diff1 = fabs(Teq_1[0]-log10(temp_elc_av));
      diff2 = fabs(Teq_2[0]-log10(temp_elc_av));
      if (diff1 > diff2) t_idx += 1;
    }
 
    q_idx = m0_idx*up->resTe + t_idx;
    double *iz_dat_d = (double*) gkyl_array_fetch(ioniz_data, q_idx);
    //printf("m0 %g Te %g iz_dat %g\n", m0_elc_av, temp_elc_av, iz_dat_d[0]);
    coef_iz_d[0] = iz_dat_d[0]/cell_av_fac; 
  }
}

__global__ static void
gkyl_iz_react_rate_neut_cu_ker(const struct gkyl_dg_iz *up, const struct gkyl_range conf_rng,
  const struct gkyl_dg_prim_vars_type *calc_prim_vars_elc_vtSq,
  const struct gkyl_array *moms_elc,
  struct gkyl_array *vtSq_elc, struct gkyl_array *prim_vars_neut, struct gkyl_array *coef_iz, 
  struct gkyl_array *M0q, struct gkyl_array *Teq, struct gkyl_array *ioniz_data)
{
  int cidx[GKYL_MAX_CDIM];
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < conf_rng.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&conf_rng, tid, cidx);
    long loc = gkyl_range_idx(&conf_rng, cidx);

    const double *moms_elc_d = (const double*) gkyl_array_cfetch(moms_elc, loc);
    const double *m0_elc_d = &moms_elc_d[0];

    double *vtSq_elc_d = (double*) gkyl_array_fetch(vtSq_elc, loc);
    double *coef_iz_d = (double*) gkyl_array_fetch(coef_iz, loc); 
    calc_prim_vars_elc_vtSq->kernel(calc_prim_vars_elc_vtSq, cidx, moms_elc_d, vtSq_elc_d);
    //Find nearest neighbor for n, Te in ADAS interpolated data
    double cell_av_fac = pow(1.0/sqrt(2.0),up->cdim);
    double m0_elc_av = m0_elc_d[0]*cell_av_fac;
    double temp_elc_av = vtSq_elc_d[0]*cell_av_fac*up->mass_elc/up->elem_charge;
    double diff1 = 0;
    double diff2 = 0;
    int m0_idx, t_idx, q_idx;
    int qtotal = up->resM0*up->resTe;
    
    if (log10(m0_elc_av) < up->minLogM0) m0_idx=0;
    else if (log10(m0_elc_av) > up->maxLogM0) m0_idx=qtotal-1;
    else {
      m0_idx = (log10(m0_elc_av) - up->minLogM0)/(up->dlogM0);
      double *M0q_1 = (double*) gkyl_array_fetch(M0q, m0_idx*up->resTe);
      double *M0q_2 = (double*) gkyl_array_fetch(M0q, (1+m0_idx)*up->resTe);
      diff1 = fabs(M0q_1[0]-log10(m0_elc_av));
      diff2 = fabs(M0q_2[0]-log10(m0_elc_av));
      if (diff1 > diff2) m0_idx += 1;
    }

    if (log10(temp_elc_av) < up->minLogTe) t_idx=0;
    else if (log10(temp_elc_av) > up->maxLogTe) t_idx=up->resTe-1;
    else {
      t_idx = (log10(temp_elc_av) - up->minLogTe)/(up->dlogTe);
      double *Teq_1 = (double*) gkyl_array_fetch(Teq, t_idx);
      double *Teq_2 = (double*) gkyl_array_fetch(Teq, t_idx+1); 
      diff1 = fabs(Teq_1[0]-log10(temp_elc_av));
      diff2 = fabs(Teq_2[0]-log10(temp_elc_av));
      if (diff1 > diff2) t_idx += 1;
    }
 
    q_idx = m0_idx*up->resTe + t_idx;
    double *iz_dat_d = (double*) gkyl_array_fetch(ioniz_data, q_idx);
    //printf("m0 %g Te %g iz_dat %g\n", m0_elc_av, temp_elc_av, iz_dat_d[0]);
    coef_iz_d[0] = iz_dat_d[0]/cell_av_fac; 
  }
}

void gkyl_dg_iz_coll_elc_cu(const struct gkyl_dg_iz *up,
  const struct gkyl_array *moms_elc, const struct gkyl_array *moms_neut,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, const struct gkyl_array *b_i,
  const struct gkyl_array *f_self, struct gkyl_array *coll_iz, struct gkyl_array *cflrate)
{
  // Set auxiliary variable (b_i) for computation of gk neut prim vars
  gkyl_dg_prim_vars_transform_vlasov_gk_set_auxfields(up->calc_prim_vars_neut_gk, 
    (struct gkyl_dg_prim_vars_auxfields) {.b_i = b_i});

  gkyl_iz_react_rate_cu_ker<<<up->conf_rng->nblocks, up->conf_rng->nthreads>>>(up->on_dev, *up->conf_rng,
    up->calc_prim_vars_elc_vtSq->on_dev, up->calc_prim_vars_neut_gk->on_dev, 
    moms_elc->on_dev, moms_neut->on_dev,
    up->vtSq_elc->on_dev, up->prim_vars_neut->on_dev, up->coef_iz->on_dev,
    up->M0q->on_dev, up->Teq->on_dev, up->ioniz_data->on_dev);

  // Calculate vt_sq_iz 
  gkyl_array_copy_range(up->vtSq_iz, up->vtSq_elc, *up->conf_rng);
  gkyl_array_scale_range(up->vtSq_iz, 1/2.0, *up->conf_rng);
  gkyl_array_shiftc(up->vtSq_iz, -up->E*up->elem_charge/(3*up->mass_elc)*pow(sqrt(2),up->cdim), 0);

  // Set fmax moments
  gkyl_array_set_offset_range(up->prim_vars_fmax, 1., up->prim_vars_neut, 0, *up->conf_rng);
  gkyl_array_set_offset_range(up->prim_vars_fmax, 1., up->vtSq_iz, up->cbasis->num_basis, *up->conf_rng);

  // POM
  gkyl_proj_gkmaxwellian_on_basis_prim_mom(up->proj_max, up->phase_rng, up->conf_rng, moms_elc,
    up->prim_vars_fmax, bmag, jacob_tot, up->mass_elc, up->fmax_iz);

  // copy, scale and accumulate
  gkyl_array_set_range(coll_iz, 2.0, up->fmax_iz, *up->phase_rng);
  gkyl_array_accumulate_range(coll_iz, -1.0, f_self, *up->phase_rng);
  
  // weak multiply
  gkyl_dg_mul_op_range(*up->cbasis, 0, up->coef_iz, 0, up->coef_iz, 0, moms_neut, up->conf_rng);
  // coll_iz = coef_iz*coll_iz
  gkyl_dg_mul_conf_phase_op_range(up->cbasis, up->pbasis, coll_iz, up->coef_iz, coll_iz,
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

void gkyl_dg_iz_coll_ion_cu(const struct gkyl_dg_iz *up,
  const struct gkyl_array *moms_elc, const struct gkyl_array *moms_neut,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, const struct gkyl_array *b_i,
  struct gkyl_array *coll_iz, struct gkyl_array *cflrate)
{
  // Set auxiliary variable (b_i) for computation of gk neut prim vars
  gkyl_dg_prim_vars_transform_vlasov_gk_set_auxfields(up->calc_prim_vars_neut_gk, 
    (struct gkyl_dg_prim_vars_auxfields) {.b_i = b_i});

  gkyl_iz_react_rate_cu_ker<<<up->conf_rng->nblocks, up->conf_rng->nthreads>>>(up->on_dev, *up->conf_rng,
    up->calc_prim_vars_elc_vtSq->on_dev, up->calc_prim_vars_neut_gk->on_dev, 
    moms_elc->on_dev, moms_neut->on_dev,
    up->vtSq_elc->on_dev, up->prim_vars_neut->on_dev, up->coef_iz->on_dev,
    up->M0q->on_dev, up->Teq->on_dev, up->ioniz_data->on_dev);

  gkyl_proj_gkmaxwellian_on_basis_prim_mom(up->proj_max, up->phase_rng, up->conf_rng, moms_elc,
    up->prim_vars_neut, bmag, jacob_tot, up->mass_elc, coll_iz);

   // coef_iz = ne*<v_e*sigma>
  gkyl_dg_mul_op_range(*up->cbasis, 0, up->coef_iz, 0, up->coef_iz, 0, moms_neut, up->conf_rng);
  // coll_iz = coef_iz*coll_iz
  gkyl_dg_mul_conf_phase_op_range(up->cbasis, up->pbasis, coll_iz, up->coef_iz, coll_iz,
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

void gkyl_dg_iz_coll_neut_cu(const struct gkyl_dg_iz *up,
  const struct gkyl_array *moms_elc, const struct gkyl_array *moms_neut,
  const struct gkyl_array *f_self,
  struct gkyl_array *coll_iz, struct gkyl_array *cflrate)
{

  gkyl_iz_react_rate_neut_cu_ker<<<up->conf_rng->nblocks, up->conf_rng->nthreads>>>(up->on_dev, *up->conf_rng,
    up->calc_prim_vars_elc_vtSq->on_dev, 
    moms_elc->on_dev,
    up->vtSq_elc->on_dev, up->prim_vars_neut->on_dev, up->coef_iz->on_dev,
    up->M0q->on_dev, up->Teq->on_dev, up->ioniz_data->on_dev);

  gkyl_array_set_range(coll_iz, -1.0, f_self, *up->phase_rng);
    
  // coef_iz = ne*<v_e*sigma>
  gkyl_dg_mul_op_range(*up->cbasis, 0, up->coef_iz, 0, up->coef_iz, 0, moms_neut, up->conf_rng);
  // coll_iz = coef_iz*coll_iz
  gkyl_dg_mul_conf_phase_op_range(up->cbasis, up->pbasis, coll_iz, up->coef_iz, coll_iz,
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

struct gkyl_dg_iz*
gkyl_dg_iz_cu_dev_new(struct gkyl_rect_grid* grid, struct gkyl_basis* cbasis, struct gkyl_basis* pbasis,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng, 
  double elem_charge, double mass_elc, enum gkyl_dg_iz_type type_ion)
{
  gkyl_dg_iz *up = (struct gkyl_dg_iz*) gkyl_malloc(sizeof(struct gkyl_dg_iz));

  int cdim = cbasis->ndim;
  int pdim = pbasis->ndim;
  int poly_order = cbasis->poly_order;
  up->cbasis = cbasis;
  up->pbasis = pbasis;
  up->cdim = cdim;
  up->use_gpu = true;
  up->conf_rng = conf_rng;
  up->phase_rng = phase_rng; 

  up->elem_charge = elem_charge;
  up->mass_elc = mass_elc;

  int resM0=50, resTe=100, qpoints=resM0*resTe;
  double minM0 = 1e15, maxM0 = 1e20;
  double minTe = 1.0, maxTe = 4e3;

  up->minLogM0 = log10(minM0);
  up->minLogTe = log10(minTe);    
  up->maxLogM0 = log10(maxM0);
  up->maxLogTe = log10(maxTe);
  
  double dlogM0 = (up->maxLogM0 - up->minLogM0)/(resM0-1);
  double dlogTe = (up->maxLogTe - up->minLogTe)/(resTe-1);
  up->resM0=resM0;
  up->resTe=resTe;
  up->dlogM0=dlogM0;
  up->dlogTe=dlogTe;
  
  up->M0q = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, qpoints); 
  up->Teq = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, qpoints); 
  up->ioniz_data = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, qpoints);

  up->E = 13.6;

  // allocate fields for prim mom calculation
  up->prim_vars_neut = gkyl_array_cu_dev_new(GKYL_DOUBLE, cbasis->num_basis, up->conf_rng->volume);
  up->vtSq_elc = gkyl_array_cu_dev_new(GKYL_DOUBLE, cbasis->num_basis, up->conf_rng->volume);
  up->vtSq_iz = gkyl_array_cu_dev_new(GKYL_DOUBLE, cbasis->num_basis, up->conf_rng->volume);
  up->prim_vars_fmax = gkyl_array_cu_dev_new(GKYL_DOUBLE, 2*cbasis->num_basis, up->conf_rng->volume);
  up->coef_iz = gkyl_array_cu_dev_new(GKYL_DOUBLE, cbasis->num_basis, up->conf_rng->volume);
  up->fmax_iz = gkyl_array_cu_dev_new(GKYL_DOUBLE, pbasis->num_basis, up->phase_rng->volume);

  up->calc_prim_vars_neut_gk = gkyl_dg_prim_vars_transform_vlasov_gk_new(cbasis, pbasis, up->conf_rng, "prim", true);
  up->calc_prim_vars_elc_vtSq = gkyl_dg_prim_vars_gyrokinetic_new(cbasis, pbasis, "vtSq", true);
  
  up->proj_max = gkyl_proj_maxwellian_on_basis_new(grid, cbasis, pbasis, poly_order+1, true);

  // copy the host struct to device struct
  struct gkyl_dg_iz *up_cu = (struct gkyl_dg_iz*) gkyl_cu_malloc(sizeof(struct gkyl_dg_iz));
  gkyl_cu_memcpy(up_cu, up, sizeof(struct gkyl_dg_iz), GKYL_CU_MEMCPY_H2D);

  // set parent on_dev pointer
  up->on_dev = up_cu;
  
  return up;
}
